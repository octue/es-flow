function [sec] = preprocess_adcp_partrac(files_in, file_out, tag)
%PREPROCESS_ADCP_PARTRAC preprocesses ADCP data from Partrac into OAS *.adcp
%format, and returns a report section on data preprocessing.
%
% Inputs:
%
%       files_in        cell        Cell array of strings containing the
%                                   consecutive filenames of data files issued
%                                   by Partrac, which is usually in the form of
%                                   a series of .mat files.
%
%       file_out        string      Path and name to the output ADCP file.
%                                   The extension *.adcp will be added if not
%                                   present.
%
%       tag             string      Text tag of the instrument, which will be
%                                   used as a title of the subsection report
%                                   output.
%
% Example:
%
%       files_in = {'/Volumes/TiME Data/Pentland/Partrac Raw Data/P1464 - Pentland AD2CP/QC Datasets & Reports/P1464.03.04.02.04.D**v1 - Site 4 - Part *.mat'}
%       file_out = '/Volumes/TiME Data/Pentland/ADCPs/TiME Nortek AD2CP'
%       sec = preprocess_adcp_partrac(folder_in, pattern_in, file_out);
%
% Outputs:
%
%       sec             [1 x 1]     Handle to a report Section object on data
%                                   preprocessing.
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


% Add extension to the output file
if ~strcmpi(file_out(end-4:end), '.adcp')
    file_out = [file_out '.adcp'];
end

% Check we're not overwriting
assert(exist(file_out, 'file')==0, 'Output *.adcp file already exists. Terminating to avoid accidental overwrite.')

% Check the input pattern and folder
for i = 1:numel(files_in)
    assert(exist(files_in{i}, 'file')~=0, 'One or more input files does not exist.')
end

% We load and stitch the partrac data. This also adds zUnit and direction
% fields, and corrects for timebase drift
[adcpData, items] = loadAndStitchPartrac(tag, files_in{:}); %#ok<ASGLU>

% Initialise a report section
sec = SubSection(1, tag, items{:});
    
% We save the data to a -v7.3 file format
save(file_out, '-struct', 'adcpData', '-v7.3')


end

function [adcpData, items] = loadAndStitchPartrac(tag, varargin)
% Load partrac file format ADCP data and stitch multiple files together

    % Number of files to stitch
    nFiles = nargin-1;
    
    % Create an empty data structure to start off with
    adcpData.u = [];
    adcpData.v = [];
    adcpData.w = [];
    adcpData.z = [];
    adcpData.d = [];
    adcpData.t = [];
    adcpData.flags.spikes = [];
    
    % Successively load each file, overwriting s each time to minimise memory
    % usage
    nBefore = 0;
    nPoints = zeros(nFiles,1);
    for iFile = 1:nFiles
        
        % Get the name and path, display progress
        [~, name] = fileparts(varargin{iFile});
        dispnow(['Loading file ' name])
        
        % Load the current chunk of data
        s = load(varargin{iFile});
        
        % Standard variables transformed to the OAS data format
        adcpData.u = [adcpData.u s.East'];
        adcpData.v = [adcpData.v s.North'];
        adcpData.w = [adcpData.w s.Vertical'];
        if iFile == 1
            adcpData.z = [adcpData.z; s.BinHeight(:)];
        end
        
        adcpData.d = [adcpData.d s.WaterDepth(:)'];
        
        % For this dataset, we have a ManualSpikes field where John has manually
        % removed data. We'll make sure to retain it but take care when
        % truncating the dataseries as the indices won't work any more.
        spikes = false;
        if isfield(s, 'ManualSpike')
            adcpData.flags.spikes = [adcpData.flags.spikes; bsxfun(@plus,s.ManualSpike, [nBefore 0])];
            spikes = true;
        end
       
        % Add the offset used to keep the spike indices correct for the next
        % loop iteration
        nBefore = nBefore + numel(s.WaterDepth);
        
        % We also have metadata to load and store
        adcpData.flags.metadata(iFile) = s.Metadata;
        
        % Add times
        adcpData.t = [adcpData.t s.DateStamp'];
        
        % Note the number of points added by each file ready for rebasing in the
        % event that data is not monotonic.
        nPoints(iFile) = numel(s.DateStamp);
        
%         % Manual check on the non-monotonic data spike
%           if iFile == 1
%               datevec(s.DateStamp(1))
%               datevec(s.DateStamp(2))
%           end
%         if iFile == 3
%             date1 = s.DateStamp(end);
%             datevec(s.DateStamp(end))
%             
%         elseif iFile == 4
%             date2 = s.DateStamp(1);
%             datevec(s.DateStamp(1))
%             disp('here')
%             datevec(date2-date1)
%             
%         end
    end
    
    % Check for monotonic times
    monotonic = true;
    if ~isMonotonic(adcpData,double(10*eps('single')))
        OASFigures
        fh = raiseFigure(['Monotonic Data Check - ' tag]);
        clf
        dt = datevec(diff(adcpData.t));
        subplot(2,1,1)
        plot(dt(:,6))
        xlabel('Data index')
        ylabel('\Delta t (s)')
        ax1 = gca;
        subplot(2,1,2)
        plot(adcpData.t-adcpData.t(1))
        xlabel('Data index')
        ylabel('t (s)')
        ax2 = gca;
        linkaxes([ax1,ax2],'x')
        monotonic = false;
    end
    
    % Manually add zUnit
    adcpData.zUnit = 0.6; % m of ADCP head above seabed

    % Manually add principal direction, default u = eastings, v = northings direction used
    adcpData.direction = 0;

    % Force monotonic timebase and account for drift
    if isfield(s.Metadata,'TimeDriftRecovery')
        drift = s.Metadata.TimeDriftRecovery; % seconds. NB CHECK THIS IS FORWARD OR BACKWARD DRIFT - ASSUMING FORWARD
    elseif isfield(s.Metadata,'TimeDriftRecoverySec')
        drift = s.Metadata.TimeDriftRecoverySec; % seconds. NB CHECK THIS IS FORWARD OR BACKWARD DRIFT - ASSUMING FORWARD
    else
        warning('No drift rate available in metadata. Applying 0 correction for drift')
        drift = 0;
    end
    if ~isnumeric(drift)
        drift = 0;
    end
    % Note that using the 'Inf' option on correctTime is not recommended unless
    % you've checked that any variation from non-monotonic data is acceptable.
    % This is OK here because we're putting the visualisation into a report
    adcpData = correctTime(adcpData,drift,Inf);
        
    % Notes on the preprocessing method
    items = {Paragraph('Data supplied and quality-controlled by the survey partner (Partrac) was converted to OAS standard *.adcp format, with velocities in the earth-fixed reference frame.')};
    items{end+1} = Paragraph('Data preprocessing included:',...
                             '\begin{itemize}',...
                                ['\item[-] Stitch of total ' num2str(nargin) ' files.'],...
                                '\item[-] Timebase continuity check (i.e. no gaps between files).',...
                                '\item[-] Clock rate consistency check (clock rate does not vary by more than $10x$ single precision or 1\% between successive data points).',...                                
                                ['\item[-] Added reference height above seabed of $' num2str(adcpData.zUnit) 'm$.'],...                                
                                ['\item[-] Timebase drift correction of $' num2str(drift) 's$ (positive drift represents clocking too fast compared to UTC).'],...
                                '\item[-] Clock rate monotonicity correction (forced a strictly monotonic timebase increase, to double precision, for purposes of numeric processing).',...
                              '\end{itemize}');
    if spikes
        items{end+1} = Paragraph('Manual spike flags denoting problem values in the dataset identified by Partrac were added.');
    else
        items{end+1} = Paragraph('No spike flags (denoting problem values in the dataset) were identified or added.');
    end
    
    if ~monotonic
        items{end+1} = Figure(fh, 0.8, 'Check of data monotonicity. Large changes in $\Delta t$ indicate significant skips in the recording. Smaller changes may be due to numerical precision issues in the instrument timebase or missed values from the dataset.');
        items{end+1} = Paragraph('Data in the input files was found to be non-monotonic. The timebase increase between each successive datapoint is shown in figure ', Ref(items{end}), '. The correction applied to account for this is to force monotonicity. The quality control process included an inspection to ensure that a forcing approach is approporiate (i.e. that the skip in timebase between datapoints is neglibible for the purposes of analysis).');
    end
    items{end+1} = Paragraph('\FloatBarrier');
        

end
