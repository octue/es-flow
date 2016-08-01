classdef InstrumentAnalysis
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
        
        function obj = InstrumentAnalysis(resultsFolder, reportFolder, name)
            
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
        
        function obj = set.DataFiles
        end
                
    end
    
end

