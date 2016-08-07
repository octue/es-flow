classdef InstrumentAnalysis < handle
    %INSTRUMENTANALYSIS Provides an analysis framework for an instrument
    %
    % Future Improvements:
    %
    %   [1] Subclassing to specific instruments
    %
    % References:
    %
    %   [1] 
    %
    % Author:                   T. H. Clark
    % Work address:             Ocean Array Systems Ltd
    %                           Hauser Forum
    %                           3 Charles Babbage Road
    %                           Cambridge
    %                           CB3 0GT
    % Email:                    tom.clark@oceanarraysystems.com
    % Website:                  www.oceanarraysystems.com
    %
    % Copyright (c) 2016 Ocean Array Systems, All Rights Reserved.

    properties
        Tag
        UserData
    end
    
    properties (SetAccess = 'protected')
        Name                % String name of the instrument (e.g. 'LiDAR 1')
        DataFiles           % A filename or cell array of flenames which contain the input data
        ReportFolder        % Where results from this analysis are stored
        ResultsFolder       % Where reports from this analysis are stored
        InstrumentType      % Subclass to alter this
    end
    
    methods
        
        function obj = InstrumentAnalysis(resultsFolder, reportFolder, name, allowOverwrite)
            
            % Make the output directories if necessary
            if exist(resultsFolder,'dir') == 0
                mkdir(resultsFolder);
            end
            if exist(reportFolder,'dir') == 0
                mkdir(reportFolder);
            end
            
            % Check the name is a string
            assert(ischar(name), 'Input ''name'' is not a character string')
            
            % Assign the object properties
            obj.ResultsFolder = resultsFolder;
            obj.ReportFolder = reportFolder;
            obj.InstrumentType = 'generic';
                        
        end
        
        function obj = SetDataFiles(obj, type, varargin)
            %SETDATAFILES Set the data files to be used for analysis.
            %
            % obj.SetDataFiles('files',filename1,...,filenameN)
            % obj.SetDataFiles('matfiles', mf)
            
            % Input checking
            assert(nargin>2, 'Insufficient number of input arguments. You must specify type and at least one data filename or MatFiles object')
            switch lower(type)
                case 'matfiles'
                    assert(nargin ==3, 'Incorrect number of input arguments for specifying Matfile or MatFiles type data')
                    assert(isa(varargin{1},'MatFile') || isa(varargin{1},'MatFiles'), 'Input object does not correspond to given type ''MatFiles''')
                case 'files'
                    for i = 1:nargin-2
                        assert(ischar(varargin{i}), 'One or more input data file names is not a string.')
                        assert(exist(varargin{1},'file') ~= 0, 'One or more input data file names does not exist.')
                    end
                otherwise
                    error('MATLAB:InstrumentAnalysis:UnknownType','Unknown datatype specified.')
            end
            
            % Set the data file property
            if numel(varargin) == 1
                obj.DataFiles = varargin{1};
            else
                obj.DataFiles = varargin;
            end
        end
        
    end
    
end

