function [adcpInput] = filterADCP(adcpInput, filter, varargin)
%FLOWDIRECTION Returns a denoised/filtered ADCP 
%
% Methods adopted here include moving averages, Savitzky-Golay filtering and a simple mean.
%
% Directions are given in degrees anticlockwise from East.
%
% Syntax:
%
%       [adcpOutput] = filterADCP(adcpInput)
%           
%       
% Inputs:
%
%       adcpInput       structure       An OAS standard adcp data structure.
%                                       See help('loadADCP') for fields
%                                       reference.
%
%       filter          string          String denoting the filter type to be
%                                       used. Valid strings are:    
%
%                           'moving'    Default
%                                       A moving average is used to smooth the
%                                       velocity data.
%   
%                           'sgolay'    A second order Savitzky Golay filter is
%                                       used to smooth the velocity data.
%       
%                           'none'      No filtering applied to the input
%                                       signal. The output is the same as the
%                                       input. Only retained to allow consistent
%                                       options if called from other functions.
%                       
%                           'wavelet'   Wavelet based denoising is applied in
%                                       the time domain only, removing
%                                       components at timescales below the
%                                       nyquist criterion applicable to the
%                                       dataset. Nature of the filtering is
%                                       variable according to selection of the
%                                       wavelet 
%
%
% Optional Inputs (Parameter-Value pairs):
%
%       Parameter       Values
%
%       'filterSpan'    [1 x 1]         Span of the window used if
%                                       filtering with either the 'sgolay' or
%                                       'window' methods.
%                                       Selected by default to span a 5 minute
%                                       window of the results. Non default
%                                       values should be odd, non-negative
%                                       integers. Setting to 0 invokes default
%                                       behaviour.
%
%       'waveletType'   string          Default 'daub4'
%                                       Specifies the type of wavelet to
%                                       be used in cases where 'filter'
%                                       parameter is selected as 'wavelet'.
%                           'daub4'         
%
%       'plotWavelet'   [1 x 1]         Default false
%                                       Setting to true results in a plot of the
%                                       wavelet based decomposition by scale
%                                       being plotted to demonstrate the levels
%                                       retained and removed.
%
% Outputs:
%
%   	adcpOutput      structure       An OAS standard adcp data structure.
%                                       See help('loadADCP') for fields
%                                       reference.
%
% Future Improvements:      none
% Other m-files required:   none
% Subfunctions:             none
% Nested functions:         none
% MAT-files required:       none
%
% Author:           T. H. Clark
% Work address:     Hauser Forum
%                   3 Charles Babbage Road
%                   Cambridge
%                   CB3 0GT
% Email:            tom.clark@oceanarraysystems.com
% Website:          www.oceanarraysystems.com
%
% Created:          21 July 2014
% Revisions:        


%% HANDLE AND CHECK INPUTS

% Parse nondefault options into the options structure
opts.filterSpan = 0;
opts.waveletType = 'daub4';
opts.plotWavelet = false;
opts = parse_pv_pairs(opts, varargin);

% Checks on size
if (numel(size(opts.filterSpan)) ~= 2) || (any(size(opts.filterSpan) ~= [1 1]))
    error('MATLAB:InvalidInput','Size of input ''filterSpan'' parameter must be [1 x 1]')
elseif (numel(size(opts.plotWavelet)) ~= 2) || (any(size(opts.plotWavelet) ~= [1 1]))
    error('MATLAB:InvalidInput','Size of input ''plot'' parameter must be [1 x 1]')
end

% Checks on type
if ~isstruct(adcpInput)
    error('MATLAB:InvalidInput','Input ''adcpInput'' must be a data structure conforming to the OAS ADCP data structure definition (see help loadADCP)')
elseif ~ischar(filter)
    error('MATLAB:InvalidInput','Input ''filter'' parameter must be a string. Try ''none'', ''sgolay'', ''moving'' or ''wavelet''')
elseif ~ischar(opts.waveletType)
    error('MATLAB:InvalidInput','Input ''waveletType'' parameter must be a string. Try ''daub4''')
elseif ~((opts.plot == true) || (opts.true == false))
    error('MATLAB:InvalidInput','Input ''plot'' parameter must be a boolean variable set to true or false')
end

% Immediate return if no filtering options requested
if strcmpi(filter,'none')
    return
end

% Checks on value and contents
if any(~isfield(adcpInput,{'u' 'v' 'w' 'z' 't'}))
    error('MATLAB:InvalidInput','Input adcpInput is missing fields. Input must conform to OAS ADCP data structure (see help loadADCP)')
end

% Set variable defaults
if opts.filterSpan == 0 && (strcmp(opts.filter,'moving') || strcmp(opts.filter,'sgolay'))
    
    % DT between profile captures (mean allows for variable dt and missing data)
    meanDT = mean(diff(adcpInput.t));
    
    % Number of profiles in a 5 minute interval (default)
    opts.filterSpan = ceil(meanDT/datenum([0 0 0 0 5 0]));
    
end
    


%% FILTER DATA

% Apply the filter to the velocities we use before determining direction:
switch opts.filter
    case 'sgolay'
        % NB accounts for non monotonic adcpInput.t
        adcpInput.u = smooth(adcpInput.t, adcpInput.u, opts.filterSpan,'sgolay',2);
        adcpInput.v = smooth(adcpInput.t, adcpInput.v, opts.filterSpan,'sgolay',2);
        adcpInput.w = smooth(adcpInput.t, adcpInput.w, opts.filterSpan,'sgolay',2);
    case 'moving'
        % NB accounts for non monotonic adcpInput.t
        adcpInput.u = smooth(adcpInput.t, adcpInput.u, opts.filterSpan,'moving');
        adcpInput.v = smooth(adcpInput.t, adcpInput.v, opts.filterSpan,'moving');
        adcpInput.w = smooth(adcpInput.t, adcpInput.w, opts.filterSpan,'moving');
    case 'wavelet'
        % Check for non-monotonic data
        diffT = diff(adcpInput.t);
        if any(diffT ~= diffT(1))
            error('MATLAB:InvalidInput','Where wavelet based filtering methods are specified, adcp data supplied must be monotonic in time.')
        end
        
        % Apply wavelet filter
        warning('MATLAB:TODO','Wavelet functionality not implemented yet. No filtering applied.')
        
        % Plot if desired
        if opts.plotWavelet
            
        end
        
    otherwise
        error('MATLAB:InvalidInput','Input ''filter'' parameter set to an unrecognised string. Try ''none'', ''wavelet'', ''sgolay'' or ''moving''')
end






end % end main function