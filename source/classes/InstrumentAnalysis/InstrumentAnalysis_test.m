classdef InstrumentAnalysis_test < TSSTest
    %INSTRUMENTANALYSIS_TEST Test class for the InstrumentAnalysis class
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
        
        % Inherited from superclass:
        %     DiaryFile       % Full filename of the diary log file for the current set of unit tests. All command line output is echoed to this file.
        %     TestPath        % Path to the current working test directory, in which a directory structure and results files will be built.
        %     DataPath        % A path to an arbitrary common data input folder containing e.g. third party or experimental results against which individual tests can evaluate and reference.
        %     ComparePath     % A path to previously existing test harness results files, used by tests tagged 'Compare' in order to evaluate differences in output over time. 
        %     Version         % A string which either specifies the current fully-committed remote version of the TSS repo currrently in use; or equals 'uncommitted' where changes have been made to the present version that aren't committed and pushed to remote.
        %     GitInfo         % A structure containing useful info on the current state of the git repo.
        OutputPath            % A path to the directory to which results for the specific tests will be saved
        
    end
    
    % Tearup methods
    methods(TestMethodSetup)
        
        function makeResultsDir(testCase)
            %MAKERESULTSDIR Build test result subdirectories if they don't already exist
            
            % Define the output path for test-specific results and 
            testCase.OutputPath = fullfile(testCase.TestPath,'tg-flow','analyse', 'InstrumentAnalysis');
            
            % Make the subdirectories if they're not already made
            if exist(testCase.OutputPath,'dir') == 0
                testCase.fatalAssertTrue(mkdir(testCase.OutputPath), ['Unable to create test output directory: ' testCase.OutputPath])
            end
            
        end
    end
    
    
    % Teardown methods
    methods(TestMethodTeardown)
        
    end
 
    
    % Methods to test basic unit functionality
    methods (Test, TestTags = {'Unit'})
        
        function testConstructor(testCase)
            %TESTCONSTRUCTOR Tests straightforward construction of an
            %InstrumentAnalysis object
            
            resultsDir = fullfile(testCase.OutputPath,'testConstructor_results');
            reportDir = fullfile(testCase.OutputPath,'testConstructor_report');
            ins = InstrumentAnalysis(resultsDir, reportDir, 'testConstructor');
            
        end
        
    end
    
    
    % Methods to compare outputs (quantitatively) against previous outputs
    methods (Test, TestTags = {'Compare'})
        
    end
    
    
    % Methods to compare output (qualitatively) against external data
    methods (Test, TestTags = {'Validate'})
        
    end
    
end

