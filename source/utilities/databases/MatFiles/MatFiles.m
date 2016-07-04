classdef MatFiles < dynamicprops
    %MATFILES 
    
    properties (SetAccess = protected)
        
        % Name (string) or names (cell array of strings) containing the .MAT
        % files to read
        Files
        
        % When the same field is read from multiple files, fields are
        % concatenated together along a prespecified dimension:
        CatDimension
        
    end
    
    properties (Access = private, Hidden = true)
        
        % Mat file mapping objects
        MatObj
        
        % Information structure of size [1 x NFields] containing the following fields:
        %   .Name       Name of field as stored in .MAT files.
        %   .Alias      Name of field as accessed from the MatFiles object.
        %               Defaults to the same as .Name.
        %   .Size       Dimensions of the array which would result if this field
        %               were taken from each file, and concatenated along
        %               CatDimensions from all the individual fields in each
        %               file.
        %   .Class      Class of the array stored in each field. Note that each
        %               field must have the same class in all files.
        %   .Bytes      The size in bytes of and array having .Size size and
        %               .Class numeric type (i.e. the memory required to hold
        %               the concatenated field from all the MAT files).
        %   .fileInds   [nFiles x 1] Cumulative sum of entries in each file used
        %               to determine which parts of which files to read when
        %               indexing into a 'Big Data' array.
        FieldInfo
        
        % Allowed fields (cell array of strings). This allows us to create a
        % MatFiles database with only a subset of the fields within the .MAT
        % files.
        AllowedFields
        
        % Structure for retaining the information about the files added to the
        % database
        FieldStruct = struct('name',[],'size',[],'bytes',[],'class',[],'global',[],'sparse',[],'complex',[],'nesting',[],'persistent',[])
        
        % Cell array to contain the meta.dynamicproperty handle objects
        MdpHandles
        
    end
    
    methods
        
        function mfo = MatFiles(fileList, varargin)
            
            % If field names are parsed, then use the specified list, otherwise
            % adopt fieldnames from the first .MAT file.
            if nargin > 1
                
                % Check validity and add the field names to the object
                if ischar(varargin{1})
                    mfo.AllowedFields{1} = varargin{1};
                else
                    mfo.AllowedFields = varargin{1};
                end
                if iscell(mfo.AllowedFields)
                    for iField = 1:numel(mfo.AllowedFields)
                        if ~ischar(mfo.AllowedFields{iField})
                            error('MatFiles:InvalidInput','allowedFields must be a string or a cell array of strings')
                        end
                    end
                else
                    error('MatFiles:InvalidInput','allowedFields must be a string or a cell array of strings')
                end
                
            end
            
            % Add all the files in the list to the database
            if ~iscell(fileList)
                error('MatFiles:InvalidInput','fileList must be a cell array of strings. If you want to use a single file, put it in a 1x1 cell; or try matlab''s own matfile() instead')
            end
            for iFile = 1:numel(fileList)
                mfo = AddFile(mfo, fileList{iFile});
            end
            
            % Add the custom or default concatenation directions
            if nargin > 2
                mfo.CatDimension = varargin{2};
            else
                mfo.CatDimension = ones(1,numel(mfo.FieldInfo));
            end
            
            % Update the FieldInfo from the FieldStruct and add dynamic
            % properties
            mfo = SetFields(mfo);
            
        end 

        function flds = fieldnames(mfo, varargin)
            %FIELDNAMES Get structure field names or object properties.
            %   NAMES = FIELDNAMES(S) returns a cell array of strings containing 
            %   the names of the fields in structure S.
            %
            %   NAMES = FIELDNAMES(Obj) returns a cell array of strings
            %   containing the names of the public properties of Obj. MATLAB
            %   objects can overload FIELDNAMES and define their own behavior.
            %   Where Obj is of class MatFiles, NAMES correspond to the aliased
            %   names in Obj, not the original fieldnames.
            %
            %   NAMES = FIELDNAMES(Obj,'-mat') where Obj is a MatFiles object
            %   returns a cell array of strings corresponding to the actual
            %   field names stored in the database of files; rather than their
            %   aliases (as returned by NAMES = FIELDNAMES(Obj).

            if (nargin > 1) && strcmpi(varargin{1},'-mat')
                % Return field names as stored in the mat files.
                flds = {mfo.FieldInfo(:).Name};
            else
                % Return field names as 'seen' by the user (i.e. the aliases,
                % not the mat file contents)
                flds = {mfo.FieldInfo(:).Alias};   
            end
        end
        
        function data = ReadData(mfo, fieldIdx, S) 
            
            % Initialise the output as empty
            data = [];
            
            % The concatenating dimension is:
            catDim = mfo.CatDimension(fieldIdx);
            
            % Get the indices of the files that contain each requested element 
            fileInds = mfo.FieldInfo(fieldIdx).FileInds;
            
            % Get the size of the output variable
            requestedInds = S.subs{catDim};
            
            % Check for errors
            if any(requestedInds > fileInds(end)) || any(requestedInds <= 0)
                error('MatFiles:OutOfBounds','Requested index out of bounds of array')
            elseif any(requestedInds ~= round(requestedInds))
                error('MatFiles:NonintegerIndex','Requested index must be an integer value')
            end
            
            % Check other dimensions
            fieldSize = mfo.FieldInfo(fieldIdx).Size;
            if numel(S.subs) ~= numel(fieldSize)
                error('MatFiles:OutOfBounds','Requested index out of bounds of array')
            end                
            for iCheck = 1:numel(S.subs)
                if (iCheck ~= catDim) && (S.subs{iCheck} ~= fieldSize(iCheck))
                    error('MatFiles:OutOfBounds','Requested index out of bounds of array')
                end
            end
                    
                
            
            % Use fileInds to get a logical map showing which file each
            % requested index is in. NB this could get memory hungry for large
            % variables and lots of files... could be a better way...
            fileUpperInds = fileInds;
%             fileLowerInds = fileInds - fileInds(1);
            fileLowerInds = ones(size(fileInds));
            fileLowerInds(2:end) = fileLowerInds(2:end) + fileInds(1:end-1);
            
%             fileMap = bsxfun(@gt, requestedInds, fileLowerInds) & bsxfun(@le, requestedInds, fileUpperInds);
            fileMap = bsxfun(@ge, requestedInds, fileLowerInds) & bsxfun(@le, requestedInds, fileUpperInds);
            
            % Logical mask for the files to use
            fileUsed = any(fileMap,2);
            
            for fileCtr = 1:numel(fileInds)
                
                
                if fileUsed(fileCtr)
                    
                    % Get the subindices required from just this file
%                     inds = requestedInds(fileMap(fileCtr,:)) - fileLowerInds(fileCtr);
                    inds = requestedInds(fileMap(fileCtr,:)) - fileLowerInds(fileCtr) + 1;
                    fileSubs = S.subs;
                    fileSubs{catDim} = inds;
                    
                    % Read the section of data from the file
                    matObj = mfo.MatObj{fileCtr};
                    try
                    dataPart = matObj.(mfo.FieldInfo(fieldIdx).Name)(fileSubs{:});
                    catch me
                        fileUpperInds
                        fileLowerInds
                        
                        fileLowerIndsOld = fileInds - fileInds(1)
                        
                        if fileCtr>1
                            fulM1 = fileUpperInds(fileCtr-1)
                            fliM1 = fileLowerInds(fileCtr-1)
                        end
                            fUI = fileUpperInds(fileCtr)
                            fLI = fileLowerInds(fileCtr)
                        if fileCtr<numel(fileInds)
                            fulP1 = fileUpperInds(fileCtr+1)
                            fliP1 = fileLowerInds(fileCtr+1)
                        end
                    
                        fieldIdx
                        nm = mfo.FieldInfo(fieldIdx).Name
                        bytes = mfo.FieldInfo(fieldIdx).Bytes
                        sz = mfo.FieldInfo(fieldIdx).Size
                        fileCtr
                        fileInds = mfo.FieldInfo(fieldIdx).fileInds
                        filename = mfo.Files{fileCtr}
                        fsSize = size(fileSubs)
                        fs1Size = size(fileSubs{1})
                        [maxval maxind] = max(fileSubs{1})
                        [minval minind] = min(fileSubs{1})
                        %filesubstart = fileSubs{1}
                        %filesubstart = fileSubs{2}
                        %fileSubEnd = fileSubs{end}
                        disp(me)
                        error('here')
                    end
                    
                    % This could be preallocated rather than grown... which
                    % would be much faster for big variables or where there are
                    % lots of concatenations but for now let's just get the code
                    % completed and working.
                    data = cat(catDim, data, dataPart);
                    
                end
            end
                
            
        end
        
        function mfo = SetAlias(mfo, fieldName, alias)
            %SETALIAS Renames a field using an alias. The original field name is
            % retained and used internally to access that field in the .MAT
            % files (i.e. the data in the MAT files aren't renamed) but are
            % simply accessed differently.
            %
            % Inputs fieldName and alias can be single strings or cell arrays of
            % strings. If cell arrays, they must be the same size; fields having
            % name fieldName{i} are renamed to alias{i}.
            
            if ~strcmp(class(fieldName),class(alias)) || ~(ischar(fieldName) || iscell(fieldName))
                error('MatFiles:InvalidInput','Input fieldName and alias must be of the same type, and either strings or cell arrays of strings')
            end
            
            if ischar(fieldName)
                fieldName = {fieldName};
                alias = {alias};
            end
            
            % For each of the input pairs
            for i = 1:numel(fieldName)
                
                % Get the field index of this name (throws an error if there's
                % no matching name)
                fieldIdx = GetFieldIndexFromName(mfo, fieldName{i});
                
                % Check that the new alias isn't the same as another alias
                % presently being used (unless we're overwriting the alias with
                % itself which doesn't do anything but shouldn't throw an error)
                currentAliases = fieldnames(mfo);
                currentFields  = fieldnames(mfo,'-mat');
                currentAliases(fieldIdx) = [];
                currentFields(fieldIdx)  = [];
                if any(ismember(currentAliases,alias{i}))
                    error('MatFiles:ConflictingAlias','Attempt to set an alias already in use as a field name')
                elseif any(ismember(currentFields,alias{i}))
                    warning('MatFiles:AliasConflictsWithFieldName',['Alias ' alias{i} ' conflicts with name of a different field in the database files. Aliases will be used to reference the data but this mapping may lead to confusion'])
                end
                
                % Set the field info to contain the appropriate alias
                mfo.FieldInfo(fieldIdx).Alias = alias{i};
                
                % Set the name of the dynamic property to the alias by tearing
                % down then instantiating new dummy variable.
                delete(mfo.MdpHandles{fieldIdx});
                mfo = AddProperty(mfo, fieldIdx);
                
            end
            
        end
        
        function mfo = ClearAliases(mfo)
            %CLEARALIASES Clears all aliases and resets field names.
            fldNames = fieldNames(mfo,'-mat');
            mfo = SetAlias(mfo,fldNames,fldNames);
        end
        
    end % End public methods
        
    methods (Hidden = true)
        
        function mfo = AddProperty(mfo, fieldIdx)
            %ADDFIELD Adds dynamic properties allowing fields in the database to
            % be accessed just like fields of a structure.
            
            % Get the alias, which will be used as the property name
            alias = GetAlias(mfo, fieldIdx);
            mdp = addprop(mfo, alias);
            mdp.SetAccess = 'protected';
            
            % Create a dummy variable to hold size and class details. Also can
            % be referenced and indexed just like a numeric datatype, as it
            % defines its own subsref.
            sz    = GetSize(mfo,  fieldIdx);
            cls   = GetClass(mfo, fieldIdx);
            bytes = GetBytes(mfo, fieldIdx);
            dv = DiskVar(mfo, fieldIdx, sz, cls, bytes);
            
            % Set the newly created property to contain this dummy variable
            mfo.(alias) = dv;
            
            % Store the mdp handle so we can alias and tear down and change
            % permissions etc later on
            mfo.MdpHandles{fieldIdx} = mdp;
            
        end
        
        function mfo = AddFile(mfo, fileName)
            %ADDFILE Adds a file to the database. mfo is the object; fileName
            %is a string containing path and file name to be added
            disp(['MatFiles: Adding file ' fileName])
            % Determine fields in the current file (whos also checks file is
            % present and a valid mat file).
            fldInfo = whos('-file', fileName);

            % Add the file name and associated MatFile Object to the database
            newInd = numel(mfo.Files) + 1;
            mfo.Files{newInd} = fileName;
            mfo.MatObj{newInd} = matfile(fileName);
            
            % We need to store the information for each field, in each file. We
            % store an [nFiles x nFields] structure with the .name .size .class
            % and .bytes fields. We can then turn the contents into useful
            % matrices of indices.
            
            % So for each field in the new file
            fldNames = {fldInfo.name};
            
            for i = 1:numel(fldNames)
                
                % Loop to test whether the current new file field matches a
                % field in the existing database
                j = 1;
                unassigned = true;
                while unassigned
                    
                    if j > numel(mfo.FieldInfo)
                        
                        % No match. If the field we're trying to add is a member
                        % of the allowed fields list, or there is no restriction
                        % on the fields that can be added, we add the new field
                        % into the database, extending FieldStruct by 1 column.
                        if isempty(mfo.AllowedFields) || any(ismember(mfo.AllowedFields,fldNames(i)))
                            mfo.FieldStruct(newInd,j) = fldInfo(i);
                            % Also add entry to the FieldInfo structure. Aliases
                            % are applied later, here set to their default
                            mfo.FieldInfo(j).Name  = fldInfo(i).name;
                            mfo.FieldInfo(j).Alias = fldInfo(i).name;
                        end
                        
                        % Either way, we're done with this new field.
                        unassigned = false;
                            
                    elseif strcmp(fldNames{i},mfo.FieldInfo(j).Name)
                            
                        % Match. Add the field from the incoming file into
                        % the correct column of the already-existing FieldStruct.
                        mfo.FieldStruct(newInd,j) = fldInfo(i);
                        unassigned = false;
                            
                        % TODO It might be easier to directly create the
                        % fileInds and accumulate the sizes and do checks etc
                        % here than in the SetFields() method
                    else
                        
                        % No match found yet. Increment to check the next field
                        % in the existing database
                        j = j+1;
                    end
                    
                end
                
            end
            
        end
        
        function mfo = SetFields(mfo)
            %SETFIELDS determines field sizes and sets up dummy variables
            
            % Update FieldInfo from the current FieldStruct
            for iField = 1:size(mfo.FieldStruct,2)
                
                % Get the number of entries along the concatenating direction in
                % each file; then accumulate this so we can determine which
                % files to look in to read data of a particular index
                fileInds = zeros(size(mfo.FieldStruct,1),1);
                for iFile = 1:size(mfo.FieldStruct,1)
                    fileInds(iFile,1) = mfo.FieldStruct(iFile,iField).size(mfo.CatDimension(iField));
                end
                fileInds = cumsum(fileInds);
                
                % Get the fields and aliases which we use for reference if we
                % throw errors                
                flds = fieldnames(mfo,'-mat');
                aliases = fieldnames(mfo);
                
                % Get the size of the stitched together 'Big Data' array
                sizes = cell2mat({mfo.FieldStruct(:,iField).size}');
                sz = zeros(1,size(sizes,2));
                for iDim = 1:size(sizes,2)
                    if (iDim ~= mfo.CatDimension(iField))
                        
                        % Gives the size of the stitched 'Big Data' array in the
                        % current dimension
                        sz(1,iDim) = sizes(1,iDim);
                        
                        % Check: If sizes vary down the column then the
                        % concatenation won't work and we should pick this up
                        % BEFORE reading the data from disc and trying to
                        % concatenate it!
                        if any(sizes(:,iDim) ~= sz(1,iDim))
                            disp('Sizes of fields to be concatenated:')
                            disp(sizes)
                            error(['Dimensions along which you are not concatenating are not consistent for MAT file field ' flds{iField} ' (alias ' aliases{iField} ')'])
                        end
                    else
                        % The size of the stitched 'Big Data' array in the
                        % concatenating dimension is the sum of all the sizes
                        % along this dimension in all the files (from above)
                        sz(1,iDim) = fileInds(end);
                    end
                end
                
                % Get the size in memory that the stitched together big data
                % array would occupy
                bytes = sum(cell2mat({mfo.FieldStruct(:,iField).bytes}));
                
                % Get the classes in all the fields and throw an error if there
                % are mismatched numeric classes
                classCell = {mfo.FieldStruct(:,iField).class}';
                for i = 1:numel(classCell)
                    if ~isempty(classCell{i}) && ~strcmp(classCell{i},classCell{1})
                        disp(['Field Class: ' classCell{1}])
                        disp(['Class in file ' mfo.Files{i} ':  ' classCell{i}])
                        error(['Mismatched data class between different files for field "' flds{iField} '" (alias "' aliases{iField} '")'])
                    end
                end
                class = classCell{1};
                
                % Update the FieldInfo for this field. NB aliases and names
                % don't get updated - they're already there (used above) and
                % updates are handled elsewhere.
                mfo.FieldInfo(iField).Size = sz;
                mfo.FieldInfo(iField).Bytes = bytes;
                mfo.FieldInfo(iField).Class = class;
                mfo.FieldInfo(iField).FileInds = fileInds;
                
            end
            
            % Build up the dynamic properties with dummy variables
            for fieldIdx = 1:numel(mfo.FieldInfo)
                mfo = AddProperty(mfo, fieldIdx);
            end
            
        end
        
        function sz = GetSize(mfo, fieldIdx)
            sz = mfo.FieldInfo(fieldIdx).Size;
        end
        
        function cls = GetClass(mfo, fieldIdx)
            cls = mfo.FieldInfo(fieldIdx).Class;
        end
        
        function bytes = GetBytes(mfo, fieldIdx)
            bytes = mfo.FieldInfo(fieldIdx).Bytes;
        end
        
        function alias = GetAlias(mfo, fieldIdx)
            alias = mfo.FieldInfo(fieldIdx).Alias;
        end
        
        function fn = GetFieldName(mfo, fieldIdx)
            fn = mfo.FieldInfo(fieldIdx).Name;
        end
        
        function fieldIdx = GetFieldIndexFromAlias(mfo, alias)
            % GETFIELDINDEXFROMALIAS Find index of a field by alias
            
            % Get alias list
            fldAliases = fieldnames(mfo);
            
            % Use ismember to get a logical falue
            fieldIdx = find(ismember(fldAliases, alias) == true);
            
            % If we reach here then there's no matching alias
            if isempty(fieldIdx)
                error('MatFiles:InvalidAlias','Alias does not correspond to any field in the current MatFiles Object')
            end
           
        end
        
        function fieldIdx = GetFieldIndexFromName(mfo, fieldName)
            % GETFIELDINDEXFROMALIAS Find index of a field by mat file fieldname
            
            % Get alias list
            flds = fieldnames(mfo,'-mat');
            
            % Find its index
            fieldIdx = find(ismember(flds, fieldName) == true);
            
            % If we reach here then there's no matching name
            if isempty(fieldIdx)
                error('MatFiles:InvalidName','Field name does not correspond to any field in the current MatFiles Object')
            end
           
        end
        
    end % End hidden methods
    
end

