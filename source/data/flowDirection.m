function [direction] = flowDirection(adcpInput, varargin)
%FLOWDIRECTION Determines direction of motion for flow in ADCP data
% Not as easy as it sounds. Each bin may have a different direction within the
% same profile, and direction changes between successive profiles in time. The
% direction is influenced by noise, waves and turbulence.
%
% To resolve the typical direction throughout the water column, it is possible
% to use mean direction over all the bins, or correct the direction for volume
% flux (i.e. bins with higher velocity contribute more strongly to the
% computation of direction). Alternatively a fixed bin may be used for the
% reference.
% 
% To resove a typical direction over time for a given bin (or the average of
% bins in the water column) some kind of filtering is required. Methods adopted
% here include moving averages, Savitzky-Golay filtering and a simple mean.
%
% Directions are given in degrees anticlockwise from East.
%
% Syntax:
%
%       [direction] = flowDirection(adcpInput)
%       [direction] = flowDirection(adcpInput, 'Parameter', value, ...)
%           
%       
% Inputs:
%
%       adcpInput       structure       An OAS standard adcp data structure.
%                                       See help('loadADCP') for fields
%                                       reference.
%
% Optional Inputs (Parameter-Value pairs):
%       
%       Parameter       Values
%
%       'filter'        'none'          Default
%                                       No filtering applied to the input
%                                       signal. Note that ascertaining 
%   
%                       'sgolay'        A second order Savitzky Golay filter is
%                                       used to smooth the adcp data before cyle
%                                       start, end and peak times are
%                                       established
%       
%                       'moving'        A moving average with window size of ten
%                                       minutes is used to determine the 
%
%       'filterSpan'    [1 x 1]         Span of the window used if
%                                       filtering with either the 'sgolay' or
%                                       'window' methods.
%                                       Selected by default to span a 10 minute
%                                       window of the results. Non default
%                                       values should be odd, non-negative
%                                       integers. Setting to 0 invokes default
%                                       behaviour.
%
%       'plot'          [1 x 1]         Boolean, Default: false
%                                       If true, plots are created showing the
%                                       signal and start/end points. If filter
%                                       options are turned on, the filtered
%                                       velocity signals used to determine bin
%                                       edges are also shown.
%                                       Plot range is selected by default to
%                                       extend 50% of the half cycle in either
%                                       direction, to show a complete cycle with
%                                       the extracted half cycle centred.
%
%       'binIndex'      [1 x 1]         Default 0
%                                       Index of the bin used to take the
%                                       direction of the flow. 
%                                       If set to 0 (the default) the average
%                                       direction of the flow throughout the
%                                       water column is used rather than an
%                                       individual bin.
%
% Outputs:
%
%   	direction       [1 x nT]        Direction in degrees anticlockwise
%                                       from East of the flow
%
% See Also: 
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
% Revision History:        	21 July 2014        Created
%                           06 April 2015       Modified documentation to OAS
%                                               header standard and solved bug
%                                               in the input checking logic
%
% Copyright (c) 2014-2015 Ocean Array Systems, All Rights Reserved.



%% HANDLE AND CHECK INPUTS

% Parse nondefault options into the options structure
opts.filter = 'none';
opts.filterSpan = 0;
opts.plot = false;
opts.binIndex = 0;
opts = parse_pv_pairs(opts, varargin);

% Checks on size
if (numel(size(opts.filterSpan)) ~= 2) || (any(size(opts.filterSpan) ~= [1 1]))
    error('MATLAB:InvalidInput','Size of input ''filterSpan'' parameter must be [1 x 1]')
elseif (numel(size(opts.binIndex)) ~= 2) || (any(size(opts.binIndex) ~= [1 1]))
    error('MATLAB:InvalidInput','Size of input ''binIndex'' parameter must be [1 x 1]')
elseif (numel(size(opts.plot)) ~= 2) || (any(size(opts.plot) ~= [1 1]))
    error('MATLAB:InvalidInput','Size of input ''plot'' parameter must be [1 x 1]')
end

% Checks on type
if ~isstruct(adcpInput) && ~isa(adcpInput,'matlab.io.MatFile')
    error('MATLAB:InvalidInput','Input ''adcpInput'' must be a data structure conforming to the OAS ADCP data structure definition (see help loadADCP) or a matfile equivalent')
elseif ~ischar(opts.filter)
    error('MATLAB:InvalidInput','Input ''filter'' parameter must be a string. Try ''none'', ''sgolay'' or ''moving''')
elseif ~((opts.plot == true) || (opts.plot == false))
    error('MATLAB:InvalidInput','Input ''plot'' parameter must be a boolean variable set to true or false')
end

% Checks on value and contents
if any(~isfield(adcpInput,{'u' 'v' 'w' 'z' 't'}))
    if any(~[isprop(adcpInput,'u') isprop(adcpInput,'v') isprop(adcpInput,'w') isprop(adcpInput,'z') isprop(adcpInput,'t')])
        error('MATLAB:InvalidInput','Input adcpInput is missing fields. Input must conform to OAS ADCP data structure (see help loadADCP)')
    end
elseif (opts.binIndex > numel(adcpInput.z)) || (opts.binIndex < 0)
    error('MATLAB:InvalidInput','Input ''binIndex'' parameter is out of range and should either index an entry in adcpInput.z, or be set = 0 (resulting in default behaviour of averaging all bins to determine direction of flow)')
end

% Set variable defaults
if opts.filterSpan == 0 && ~strcmp(opts.filter,'none')
    
    % DT between profile captures (mean allows for variable dt and missing data)
    meanDT = median(diff(adcpInput.t));
    
    % Number of profiles in a 10 minute interval (default)
    opts.filterSpan = ceil(fs*10*60/meanDT);
    
end


%% GET DIRECTION OF FLOW

% Default behaviour is to average down all the profiles, 
switch opts.binIndex
    case 0
        % Average down the bin profiles
        u = nanmean(adcpInput.u,1);
        v = nanmean(adcpInput.u,1);
    otherwise
        u = adcpInput(opts.binIndex,:);
        v = adcpInput(opts.binIndex,:);
        if isnan(u) || isnan(v)
            warning('NaN values exist in this bin. Choose a lower bin or use the removeNaNs function.')
        end
end

% Apply the filter to the velocities we use before determining direction:
switch opts.filter
    case 'none'
        uFilt = u;
        vFilt = v;
    case 'sgolay'
        % NB accounts for non monotonic adcpInput.t
        uFilt = smooth(adcpInput.t, u,opts.filterSpan,'sgolay',2);
        vFilt = smooth(adcpInput.t, v,opts.filterSpan,'sgolay',2);
    case 'moving'
        % NB accounts for non monotonic adcpInput.t
        uFilt = smooth(adcpInput.t, u,opts.filterSpan,'moving');
        vFilt = smooth(adcpInput.t, v,opts.filterSpan,'moving');
    otherwise
        error('MATLAB:InvalidInput','Input ''filter'' parameter set to an unrecognised string. Try ''none'', ''sgolay'' or ''moving''')
end

direction = atan2d(vFilt, uFilt);




%% GRAPHICAL OUTPUT

if opts.plot
    raiseFigure('Direction of Flow')
    clf
    plot(adcpInput.t, direction, 'b-')
    xlabel('Time [datenum form] (s)')
    ylabel('Degrees ACW from East')
    ylim([-180 180])
    hold on
end



end % end main function