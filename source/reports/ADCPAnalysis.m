classdef ADCPAnalysis < InstrumentAnalysis
    %ADCPANALYSIS Maps to data and provides methods
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
%         Inherited from superclass:
%         Tag
%         UserData
    end
    
    properties (SetAccess = 'protected')
%         Inherited from superclass:
%         Name                % String name of the instrument (e.g. 'LiDAR 1')
%         DataFiles           % A filename or cell array of filenames which contain the input data
%         ReportFolder        % Where results from this analysis are stored
%         ResultsFolder       % Where reports from this analysis are stored
%         InstrumentType      % Subclass to alter this
          Data                % Matfile object mapping to the data on disc (read only)
          Results             % Matfile object mapping to the results file (read-write)
          WindowLength        % Window duration in seconds
          WindowOverlap       % Fraction, typically 0 or 0.5, of window overlap.
          WindowInds          % 
    end
    
    properties (Dependent)
        ResultsFile
    end
    
    methods
        
        function obj = ADCPAnalysis(varargin)
            
            % Call superclass constructor
            obj@InstrumentAnalysis(varargin{:})

            % Check that the input is a .adcp file
            assert(ischar(obj.DataFiles), 'Single OAS format .adcp file must be used to create an ADCPAnalysis.')
            assert(strcmpi(obj.DataFiles(end-4:end),'.adcp'), 'Input data file must be in OAS *.adcp format.')

            % Create a matfiles object mapping the data (read only)
            obj.Data = matfile(obj.DataFiles,'Writable', false);

            % Create a matfiles object for data output (read-write)
            % TODO create utility function to handle number increments on filenames
            [~, name] = fileparts(obj.DataFiles);
            ctr = 1;
            exists = 1;
            while exists~=0
                res_file = fullfile(obj.ResultsFolder,[name '_ADCPAnalysis_' num2str(ctr) '.mat']);
                exists = exist(res_file, 'file');
                ctr = ctr + 1;
            end
            obj.Results = matfile(obj.ResultsFile,'Writable',true);

            % Add any metadata as read-only dynamic properties
            if isfield(obj.Data, 'metadata')
                m = obj.Data.metadata;
                fn = fieldnames(m);

                % For each item in the metadata, add the dynamic property, set
                % its value to the saved data value, and restrict write access
                for i = 1:numel(fn)
                    p = addprop(obj,fn{i});
                    p.SetAccess = 'protected';
                    obj.(fn{i}) = m.(fn{i});
                end
            end
            
        end
        
        
        function rf = get.ResultsFile(obj)
            %RESULTSFILE Get the results file to which analysis data will be
            %saved
            rf = obj.Results.Source;
        end
        
        
        function SetWindowing(obj, overlap, type, len)
            
            % Error check the overlap value
            assert(overlap < 1, 'Overlapping of 1 or more creates an infinite series of windows. Recommended overlap is 0 or 0.5')
            if overlap<0
                warning('Setting overlap < 0 spaces out windows. Some data will be excluded from the analysis. This may be intentional (e.g. to accelerate processing during early preview)')
            end
            
            % Error check and process the window size
            switch lower(type)
                case 'duration'
                    length = round(obj.frequency*len);
                    if length ~= obj.frequency*len
                        dispnow('Adjusted inexactly specified window length to nearest integer value.')
                    end
                case 'length'
                    assert((len==round(len)) && (len>0), 'Specified window length must be an integer value >= 1')
                    length = len;
                otherwise
                    error('Unknown ''type'' specification for window length. Try ''length'' (the number of data points in a window) or ''duration'' (the length of a window in seconds)')
            end
            
            % Set the properties
            obj.WindowLength = length;
            obj.WindowDuration = length/obj.frequency;
            obj.WindowOffset = round(length*(1-overlap));
            obj.WindowOverlap = 1-(obj.WindowOffset/length);
            
            % Display any correction to overlap value
            if obj.WindowOverlap ~= overlap
                dispnow(['Adjusted overlap fraction to give nearest integer index offset value. New overlap = ' num2str(obj.WindowOverlap)])
            end
            
            % Compute the window start and end indices
            nPoints = numel(obj.Data.t);
            windowStart = 1:obj.WindowOffset:nPoints;
            windowEnd = obj.WindowLength:obj.WindowOffset:nPoints;
            windowStart = windowStart(1:numel(windowEnd));
            obj.WindowInds = [windowStart(:)';windowEnd(:)'];
            
        end
        
        
        function [data, results] = Window(obj, index)
            %WINDOW gets a window from the ADCP data structure.
            
            assert(~isempty(obj.WindowInds),'Use the SetWindowing() method on this ADCPAnalysis object before attempting to retrieve windows by index.')
            assert((index > 0) && (index <= size(obj.WindowInds,2)), 'Window index out of range.')
            
            data = getWindow(obj.Data, 'IndexRange', obj.WindowInds(:,index)');
            if nargout > 1
                results = [];
            end
            
        end
        
        
        function Run(obj)
            %RUN Runs all available window analyses in parallel.
            
            parfor i = 1:size(obj.WindowInds,2)
                % Do individual window analyses
                obj.AnalyseWindow(i) %#ok<PFBNS>
            end
            
        end
        
        
        function AnalyseWindow(obj, i)
            %ANALYSEWINDOW analyse a single window.
            
            % Start and end indices of this window
            inds = obj.WindowInds(:, i)';
            
            % Get the basic ADCPdata in a structure
            data = obj.Window(i);
            
            % Get the depth and time averaged direction for this window
            dir = mean(flowDirection(data));
            
            % Rotate adcp frame of reference to the local direction
            data = rotateADCP(data, dir);
            
            % Get mean velocity profiles. Transpose for more efficient memory
            % access.
            uvwBar = [mean(data.u, 2) mean(data.v, 2) mean(data.w, 2)]';
            
            % Get spectra and estimate noise levels. Store transposed for more
            % efficient memory access
            % TODO rather than preallocating, vectorise spectrumADCP to make this use case more efficient.
            % TODO calculate beam separation frequency outside this routine and
            % use it as an input to the frequency range cutoff for baptiste's
            % analysis
            for iBin = 1:nBins
                % NB use unity velocity scaling for k1 which gives more compact
                % storage.
                fRange = [1 Inf];
                [k1U, psd(:,iBin), ~, sbs(iBin), N(iBin), K(iBin), fcut(iBin)] = spectrumADCP(data, iBin, 1, fRange);
            end
            
            % Fit analytical mean profile boundary layer parameters, unweighted.
            % For first guess we use the Song data (see gupta and clark)
            % [Pi, S0, deltac0, U10].
            x0 = [0.34 26.7 max(data.z(:)) max(uvwBar(1,:))];
            profile = fitMeanProfile(data.z(:), uvwBar(1,:)', [], x0, 'lewkowicz');
            
            % Store results to the main matfile. We set it up to store frequency
            % variation down the first dimension, bin variation down the second,
            % and window variation down the third. NB we can always squeeze the
            % dimensions down.
            if isempty(obj.Results.frequency)
                obj.Results.frequency(:,1,1)        = k1U./(2*pi);
                obj.Results.beamSeparation(1,:,1)   = sbs;
                obj.Results.fcsVer                  = fcsVer;
                obj.Results.tssVer                  = tssVer;
            end
            obj.Results.psd(:,:,i)              = psd;
            obj.Results.uvwBar(1:3,:,i)         = uvwBar;
            obj.Results.windowInds(:,i)         = inds(:);
            obj.Results.t(i)                    = data.t(1);
            obj.Results.dt                      = data.t(2) - data.t(1);
            obj.Results.bapt_fRange(1,:,i)      = fRange';
            obj.Results.bapt_N(1,:,i)           = N;
            obj.Results.bapt_K(1,:,i)           = K;
            obj.Results.bapt_FCut(1,:,i)        = fcut;
            obj.Results.lew_Pi(i)               = profile.Pi;
            obj.Results.lew_S(i)                = profile.S;
            obj.Results.lew_U1(i)               = profile.U1;
            obj.Results.lew_Utau(i)             = profile.Utau;
            obj.Results.lew_deltac(i)           = profile.deltac;
            obj.Results.lew_kappa(i)            = profile.kappa;
            obj.Results.lew_resnorm(i)          = profile.resnorm;
            
        end
        
%      	function ReFit(obj)
%       end
        
    end
    
end

