function [theStruct] = xmlReadGantt(varargin)
% PARSEXML Convert XML file to a MATLAB structure.

if nargin == 1
    % Then a filename has been passed
    filename = varargin{1};
    if exist(filename,'file') ~= 2
        % Then the file doesn't exist. Throw an error and return
        disp(['xmlReadGantt.m: Attempted to read file...   ' filename])
        error('MATLAB:InvalidFilename','File not found.')
    end
else
    % No filename has been passed. Use GUI prompt.
    def_path = retrieve_path('xmlgantt');
    [file, path] = uigetfile('*.xml','Select a Merlin *.xml file',def_path);
    if isnumeric(file)
        error('MATLAB:CancelledOperation','File-load process cancelled by user.')
    end
    filename = fullfile(path,file);
    update_default_path(path,'xmlgantt')
end

% Read the file
try
   tree = xmlread(filename);
catch %#ok<CTCH>
   error('MATLAB:InvalidFile','Failed to read XML file %s.',filename);
end

% Recurse over child nodes. This could run into problems 
% with very deeply nested trees.
try
   theStruct = parseChildNodes(tree);
catch %#ok<CTCH>
   error('MATLAB:InvalidTree','Unable to parse XML file %s.',filename);
end


% ----- Subfunction PARSECHILDNODES -----
function children = parseChildNodes(theNode)
% Recurse over node children.
children = [];
if theNode.hasChildNodes
   childNodes = theNode.getChildNodes;
   numChildNodes = childNodes.getLength;
   allocCell = cell(1, numChildNodes);

   children = struct(             ...
      'Name', allocCell, 'Attributes', allocCell,    ...
      'Data', allocCell, 'Children', allocCell);

    for count = 1:numChildNodes
        theChild = childNodes.item(count-1);
        children(count) = makeStructFromNode(theChild);
    end
end

% ----- Subfunction MAKESTRUCTFROMNODE -----
function nodeStruct = makeStructFromNode(theNode)
% Create structure of node info.

nodeStruct = struct(                        ...
   'Name', char(theNode.getNodeName),       ...
   'Attributes', parseAttributes(theNode),  ...
   'Data', '',                              ...
   'Children', parseChildNodes(theNode));

if any(strcmp(methods(theNode), 'getData'))
   nodeStruct.Data = char(theNode.getData); 
else
   nodeStruct.Data = '';
end

% ----- Subfunction PARSEATTRIBUTES -----
function attributes = parseAttributes(theNode)
% Create attributes structure.

attributes = [];
if theNode.hasAttributes
   theAttributes = theNode.getAttributes;
   numAttributes = theAttributes.getLength;
   allocCell = cell(1, numAttributes);
   attributes = struct('Name', allocCell, 'Value', ...
                       allocCell);

   for count = 1:numAttributes
      attrib = theAttributes.item(count-1);
      attributes(count).Name = char(attrib.getName);
      attributes(count).Value = char(attrib.getValue);
   end
end

