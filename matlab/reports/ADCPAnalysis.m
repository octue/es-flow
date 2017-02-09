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
          WindowDuration
          WindowOffset
          
    end
    
    properties (Dependent)
        ResultsFile
    end
    
    methods
        
        function obj = ADCPAnalysis(datafile, varargin)
            
            % Call superclass constructor and add the datafile
            obj@InstrumentAnalysis(varargin{:})
            obj.SetDataFiles('files', datafile);
            
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
            obj.Results = matfile(res_file,'Writable',true);
            
            % This forces creation of the mat file
            obj.Results.frequency = [];

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
        
        
        function AttachResults(obj, file, cleanup)
            %ATTACHRESULTS Attach to an existing results file. 
            % Setting cleanup = true deletes the existing results file.
            if nargin == 2
                cleanup = false;
            end
            assert(exist(file,'file')~=0, 'The results file to attach to must already exist.')
            
            % Map into the results file we want to attach
            newres = matfile(file,'Writable',true);
            
            % Check that the windowing applied is compatible if there are
            % already results in there
            if ~all(isequal(newres.windowInds, obj.WindowInds))
                error('Attempt to load results file with different window overlaps to current settings.')
            end
            
            % TODO more checks that it's got all the necessary fields
            
            % Delete the currently attached results file
            if cleanup
                delete(obj.Results)
            end
            
            % Attach the new file
            obj.Results = newres;
            
        end
        
        
        function SetWindowing(obj, overlap, type, len)
            
            % Error check the overlap value
            assert(overlap < 1, 'Overlapping of 1 or more creates an infinite series of windows. Recommended overlap is 0 or 0.5')
            if overlap<0
                warning('Setting overlap < 0 spaces out windows. Some data will be excluded from the analysis. This may be intentional (e.g. to accelerate processing during early preview)')
            end
            
            % Error check and process the window size
            dt = datevec(obj.Data.t(1,2)-obj.Data.t(1,1));
            dt = dt(6);
            switch lower(type)
                case 'duration'
                    length = round(len/dt);
                    if length ~= len/dt
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
            obj.WindowDuration = length*dt;
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
        
        
        function [data, results] = Window(obj, index, bins)
            %WINDOW gets a window from the ADCP data structure.
            
            assert(~isempty(obj.WindowInds),'Use the SetWindowing() method on this ADCPAnalysis object before attempting to retrieve windows by index.')
            assert((index > 0) && (index <= size(obj.WindowInds,2)), 'Window index out of range.')
            
            data = getWindow(obj.Data, 'IndexRange', obj.WindowInds(:,index)');
            if nargout > 1
                % By default, return window results for all bins
                if nargin == 2
                    bins = 1:numel(data.z);
                end
                % Get the full window out
                results.frequency = obj.Results.frequency;
                results.beamSeparation = obj.Results.beamSeparation;
                results.fcsVer = obj.Results.fcsVer;
                results.fs = obj.Results.fs;
                results.psd = obj.Results.psd(:, bins, index);
                results.uvwBar = obj.Results.uvwBar(1:3,bins,index);
                results.windowInds = obj.Results.windowInds(1:2,index);
                results.t = obj.Results.t(index,1);
                results.d = obj.Results.d(index,1);
                results.flowDirection = obj.Results.flowDirection(index,1);
                results.bapt_fRange = obj.Results.bapt_fRange(1:2,index);
                results.bapt_N = obj.Results.bapt_N(bins,index);
                results.bapt_K = obj.Results.bapt_K(bins,index);
                results.bapt_FCut = obj.Results.bapt_FCut(bins,index);
                results.lew_Pi = obj.Results.lew_Pi(index,1);
                results.lew_S = obj.Results.lew_S(index,1);
                results.lew_U1 = obj.Results.lew_U1(index,1);
                results.lew_Utau = obj.Results.lew_Utau(index,1);
                results.lew_deltac = obj.Results.lew_deltac(index,1);
                results.lew_kappa = obj.Results.lew_kappa(index,1);
                results.lew_resnorm = obj.Results.lew_resnorm(index,1);
                results.sm_lew_Pi = obj.Results.sm_lew_Pi(index,1);
                results.sm_lew_S = obj.Results.sm_lew_S(index,1);
                results.sm_lew_U1 = obj.Results.sm_lew_U1(index,1);
                results.sm_lew_Utau = obj.Results.sm_lew_Utau(index,1);
                results.sm_lew_deltac = obj.Results.sm_lew_deltac(index,1);
                results.FloodRate = obj.Results.FloodRate(index,1);
                results.dirSign = obj.Results.dirSign(index,1);
            end
            
        end
        
        
        function Run(obj)
            %RUN Runs all available window analyses in parallel.
            
            warning('Check for a missing beam angle and z unit - spectrumADCP will default to 52 degrees and 0 respectively!')
            
            % Get the number of windows
            nWindows = size(obj.WindowInds,2);
            
            
            % Run all window analyses. Running the last one first implicitly
            % defines sizes. Note: Tried to parfor this; but multiple threaads
            % writing into one file is unstable (there's no threadlocking on the
            % file access).
            obj.AnalyseWindow(nWindows)
            for i = 1:nWindows
                % Do individual window analyses
                dispnow(['Processing ADCP window ' num2str(i) ' of ' num2str(nWindows)])
                obj.AnalyseWindow(i)
            end
            
            % Add the Flood rate
            obj.AddFloodRate;
            
            % Smooth the fitted parameters
            obj.RobustSmooth;
            
        end
        
        
        function AnalyseWindow(obj, i)
            %ANALYSEWINDOW analyse a single window.
            
            % Start and end indices of this window
            inds = obj.WindowInds(:, i)';
            
            % Get the basic ADCPdata in a structure
            data = obj.Window(i);
            
            % Get the depth- and time- averaged direction for this window
            dir = mean(flowDirection(data));
            % TODO sort directions generically and apply a convention on sign;
            % f/ex flood is always positive U1 and Ebb is always negative. That
            % saves having to arbitrarily switch direction
            
            % Get the depth for this window
            d = nanmean(data.d);
            
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
            nBins = numel(data.z);
            for iBin = 1: nBins
                % NB use unity velocity scaling for k1 which gives more compact
                % storage.
                fRange = [1 Inf];
                [k1U, psd(:,iBin), ~, sbs(iBin), N(iBin), K(iBin), fcut(iBin)] = spectrumADCP(data, iBin, 1, fRange);
            end
            
            % Fit analytical mean profile boundary layer parameters, unweighted.
            % For first guess we use the Song data (see gupta and clark)
            % [Pi, S0, deltac0, U10].
            uBar = sqrt(sum((uvwBar(:,:).^2),1));
            x0 = [0.34 26.7 d max(uBar(1,~isnan(uBar(1,:))))];
            profile = fitMeanProfile(data.z(:), uBar', [], x0, 'lewkowicz');
            
            % Store results to the main matfile. We set it up to store frequency
            % variation down the first dimension, bin variation down the second,
            % and window variation down the third. NB we can always squeeze the
            % dimensions down.
            if i == 1
                obj.Results.frequency         = k1U./(2*pi);
                obj.Results.beamSeparation    = sbs;
                obj.Results.fcsVer            = fcsVer('-quiet');
%                 obj.Results.tssVer            = tssVer;
                obj.Results.fs                = 1/(data.t(2) - data.t(1));
            end
            obj.Results.psd(1:size(psd,1), 1:nBins,i)	= psd;
            obj.Results.uvwBar(1:3,1:nBins,i)       = uvwBar;
            obj.Results.windowInds(1:2,i)           = inds(:);
            obj.Results.t(i,1)                      = data.t(1);
            obj.Results.d(i,1)                      = d;
            obj.Results.flowDirection(i,1)          = dir;
            obj.Results.depthAveragedSpeed(i,1)     = sqrt(sum(nanmean(uvwBar(1:2,:,:),2).^2,1));
            obj.Results.bapt_fRange(1:2,i)          = fRange(:);
            obj.Results.bapt_N(1:nBins,i)           = N(:);
            obj.Results.bapt_K(1:nBins,i)           = K(:);
            obj.Results.bapt_FCut(1:nBins,i)        = fcut(:);
            obj.Results.lew_Pi(i,1)                 = profile.Pi;
            obj.Results.lew_S(i,1)                  = profile.S;
            obj.Results.lew_U1(i,1)                 = profile.U1;
            obj.Results.lew_Utau(i,1)               = profile.Utau;
            obj.Results.lew_deltac(i,1)             = profile.deltac;
            obj.Results.lew_kappa(i,1)              = profile.kappa;
            obj.Results.lew_resnorm(i,1)            = profile.resnorm;
        end
        
        
%         function [window, sec] = ReportWindows(obj, i, bins)
%             %REPORTWINDOW Returns a report section on a window or window
%             
%             % Table of window parameters
%             % Time, PCSprings, Dirn, uvwBar, d, lew_U1, Pi, S, U_tau, d, deltac, bapt_frange 
%             
%             % Plot of bapt_N, bapt_K varying with height, tagged with values at
%             % the bins of interest.
%             
%             % Plot of PSD in the bins of interest.
%             
%         end
        

        function AddFloodRate(obj,smoothfac)
            %ADDFLOODRATE Adds a smoothed Flood Rate metric to the results. If
            %positive, the tide is flooding, if negative it is ebbing. Beware
            %this may be out of phase with tdal velocities so cannot be used as
            %a discriminant for the ebb and flood tide.
            
            if nargin == 1
                % Use default smoothing factor of 1000 which was established by
                % trial and error during the TiME project
                smoothfac = 1000;
            end
            
            % Use penalised least square spline fit, with robust outlier
            % management
            sm_d = smoothn(obj.Results.d, smoothfac, 'robust');
            
            % Differentiate (simple forward first order)
            fr = diff(sm_d)./diff(obj.Results.t);
            obj.Results.FloodRate = [fr(1); fr];
            
        end
        
        
     	function RobustSmooth(obj, maxU1, smoothfac)
            %ROBUSTSMOOTH Uses robust spline based smoothing (\cite{Garcia2010})
            %to improve parameter space fit models. Uses criteria of U1 within
            %a specified max value (input maxU1) and Pi within 3* the median Pi
            %value.
            if nargin == 1
                maxU1 = 6; % Highest tidal speed in the world... ish.
            end
            if nargin <= 2
                % Use default smoothing factor of 1000 which was established by
                % trial and error during the TiME project
                smoothfac = 1000;
            end
            
            % TODO extend to logarithmic/exponential BL fits.
            
            % Deal with outliers based on plausible ranges of U1 and Pi
            U1 = obj.Results.lew_U1;
            Pi = obj.Results.lew_Pi;
            
            % Mask on U1 criteria. NB U1 stored such that it should always be
            % positive regardless of direction.
            mask_U1 = true(size(U1));
            mask_U1(U1>maxU1) = false;
            mask_U1(U1<0) = false;
            
            % Mask on Pi criteria
            rangePi = 3*nanmedian(abs(Pi(mask_U1)));
            mask_Pi = false(size(Pi));
            mask_Pi(~isnan(Pi)) = (Pi(~isnan(Pi))<rangePi) & (Pi(~isnan(Pi))>(-1*rangePi));
            
            % Flip the sign of U1 during ebb
            fd = obj.Results.flowDirection; % temporary local var to avoid reading from disc twice
            dir_sign = sign(fd - nanmean(fd));
            
            % STore the direction signum
            obj.Results.dirSign = dir_sign;

            % Combine the logical masks from different criteria
            mask = mask_Pi & mask_U1;
            
            % Derectify U1 and smooth the approximately sinusoidal signal.
            lew_U1 = obj.Results.lew_U1.*dir_sign;
            lew_U1(~mask) = NaN;
            obj.Results.('sm_lew_U1') = smoothn(lew_U1,smoothfac,'robust','MaxIter',500).*dir_sign;
            
            % Use dynamic fieldnames to cycle through the rest of the variables
            % we want to smooth.
            % TODO update as more analyses are added
            params = {'lew_Pi','lew_S','lew_Utau','lew_deltac','lew_kappa'};
            for i = 1:numel(params)
                param = obj.Results.(params{i});
                param(~mask) = NaN;
                warning('off','MATLAB:smoothn:InitialGuess')
                % TODO implement a check on the image processing toolbox,
                % because this function behaves differently according to whether
                % you have it installed or not.
                obj.Results.(['sm_' params{i}]) = smoothn(param,smoothfac,'robust','MaxIter',500);
                
            end
            
        end
        
        
        function CleanDirection(obj)
            
           % HACK - sorts the direction problem on the sentinel V
           
           
        end
    end
    
end

