function [adcpWindow] = getHalfCycle(adcpInput, t, degreeRange, varargin)
%GETHALFCYCLE Extracts a half cycle of tidal data from an adcp data structure.
% The half cycle of data (e.g. the flood or ebb portion of the tide) surrounding
% the input time is extracted into another, smaller data structure.
%
% The change of the tide from flood to ebb is defined as the point where
% direction of the flow crosses a boundary. 
%
% Syntax:
%
%       [adcpWindow] = getHalfCycle(adcpInput, t, degreeRange)
%           Extracts the tidal half cycle surrounding time t, using input
%           degreeRange to define the directions of flood and ebb tide.
%
%       [adcpWindow] = getHalfCycle(acdpInput, t, degreeRange, 'Parameter', Value, ...)
%           Allows additional processing and criteria to be applied, including
%           use of a filter to denoise the peak flow, start and end times (the
%           crossover point and peak detection can be affected by noise, wave
%           influence and turbulence). Plotting options are also available.
%
%       [adcpWindow, startTime, endTime] = getHalfCycle(...)
%           Returns the peak flow condition during the half cycle, as well as
%           the start, peak and end times.
%       
% Inputs:
%
%       adcpInput       structure       An OAS standard adcp data structure.
%                                       See help('loadADCP') for fields
%                                       reference.
%
%       t               [1 x 1]         Time (in MATLAB datenum form) around
%                                       which the half cycle is extracted. Note
%                                       that due to noise and unsteadiness in
%                                       the measurements, times close to slack
%                                       tide points can produce false results
%                                       (the instantaneous flow can cross the
%                                       boundary several times in a short period
%                                       surrounding slack tide) so a more robust
%                                       use case is to set t to a time somewhere
%                                       clearly in the middle of a flood or ebb
%                                       tide.
%
%       degreeRange     [1 x 2]         0 <= degreeRange(i) < 360
%                                       Two bearings, in Nautical degrees, used
%                                       to define change in direction of the
%                                       tide. Whenever the flow direction
%                                       passes through either of these bearings,
%                                       a change between ebb and flood is deemed
%                                       to have occurred. Nautical bearings are
%                                       referenced as 000 at due North, 090 due
%                                       East, 180 South, 270 West.
%                                       For example [090 270] may be used as a
%                                       range where the flow is predominantly
%                                       north-south. Alternatively, in a dogleg
%                                       channel the flood tide might flow north,
%                                       and the ebb tide might flow east. In
%                                       this case, a sensible degreeRange is
%                                       [045 225].
%
% Optional Inputs (Parameter-Value pairs):
%       
%       Parameter       Values
%
%       'filter'        'none'          Default
%                                       No filtering applied to the input
%                                       signal. Note that a spiky signal will
%                                       result in inaccuracy when determining
%                                       start and end times 
%                                       flow rate, 
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
%                                       Selected by default to span a 5 minute
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
%       'binIndex'     [1 x 1]          Default 0
%                                       Index of the bin used to take the
%                                       direction of the flow for the purposes
%                                       of establishing slack water points. 
%                                       If set to 0 (the default) the average
%                                       direction of the flow throughout the
%                                       water column is used rather than an
%                                       individual bin.
%
% Outputs:
%
%   	adcpWindow      structure       An OAS standard adcp data structure
%                                       extracted from a limited part of the
%                                       input containing just the half cycle of
%                                       interest.
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
% Created:          14 July 2014
% Revisions:        


%% HANDLE AND CHECK INPUTS

% Parse nondefault options into the options structure
opts.filter = 'none';
opts.filterSpan = 0;
opts.plot = false;
opts.binIndex = 0;
opts = parse_pv_pairs(opts,varargin);

% Checks on size
if (numel(size(t)) ~= 2) || (any(size(t) ~= [1 1]))
    error('MATLAB:InvalidInput','Size of input ''time'' must be [1 x 1]')
elseif (numel(size(degreeRange)) ~= 2) || (any(size(degreeRange) ~= [1 2]))
    error('MATLAB:InvalidInput','Size of input ''degreeRange'' must be [1 x 2]')
elseif (numel(size(opts.filterSpan)) ~= 2) || (any(size(opts.filterSpan) ~= [1 1]))
    error('MATLAB:InvalidInput','Size of input ''filterSpan'' parameter must be [1 x 1]')
elseif (numel(size(opts.binIndex)) ~= 2) || (any(size(opts.binIndex) ~= [1 1]))
    error('MATLAB:InvalidInput','Size of input ''binIndex'' parameter must be [1 x 1]')
elseif (numel(size(opts.plot)) ~= 2) || (any(size(opts.plot) ~= [1 1]))
    error('MATLAB:InvalidInput','Size of input ''plot'' parameter must be [1 x 1]')
end

% Checks on type
if ~isstruct(adcpInput)
    error('MATLAB:InvalidInput','Input ''adcpInput'' must be a data structure conforming to the OAS ADCP data structure definition (see help loadADCP)')
elseif ~ischar(opts.filter)
    error('MATLAB:InvalidInput','Input ''filter'' parameter must be a string. Try ''none'', ''sgolay'' or ''moving''')
elseif ~((opts.plot == true) || (opts.true == false))
    error('MATLAB:InvalidInput','Input ''plot'' parameter must be a boolean variable set to true or false')
end

% Checks on value and contents
if any(~isfield(adcpInput,{'u' 'v' 'w' 'z' 't'}))
    error('MATLAB:InvalidInput','Input adcpInput is missing fields. Input must conform to OAS ADCP data structure (see help loadADCP)')
end
timeRange = [min(adcpInput.t) max(adcpInput.t)];
if (t > timeRange(2)) || (t < timeRange(1))
	error('MATLAB:InvalidInput','Input ''t'' is out of the time range contained in the input ADCP data structure')
elseif (opts.binIndex > numel(adcpInput.z)) || (opts.binIndex < 0)
    error('MATLAB:InvalidInput','Input ''binIndex'' parameter is out of range and should either index an entry in adcpInput.z, or be set = 0 (resulting in default behaviour of averaging all bins to determine direction of flow)')
end

% Set variable defaults
if opts.filterSpan == 0 && ~strcmp(opts.filter,'none')
    
    % DT between profile captures (mean allows for variable dt and missing data)
    meanDT = mean(diff(adcpInput.t));
    
    % Number of profiles in a 5 minute interval (default)
    opts.filterSpan = ceil(meanDT/datenum([0 0 0 0 5 0]));
    
end
    


%% GET DIRECTION OF FLOW

% Use the flowDirection function. Note that we force plot to false as it's
% plotted below in cases where the user requests plot=true.
direction = flowDirection(adcpInput, 'filter', opts.filter, 'filterSpan', opts.filterSpan, 'plot', false, 'binIndex', opts.binIndex);

% u is Eastings, v is Northings. We use the nautical convention of 0 to 360
% degrees around the clock face, so need to convert from the output of the
% direction function
degreesNautical = (360/(2*pi))*direction;
degreesNautical(degreesNautical < 0) = degreesNautical(degreesNautical < 0) + 360;

% Get a mask true in one sector and false
degreeRange = sort(degreeRange);
rangeMask = false(size(adcpInput.t));
rangeMask(degreesNautical >  degreeRange(2)) = true;
rangeMask(degreesNautical <= degreeRange(1)) = true;

% Indices into the time vector immediately preceding a direction shift
switchInds = find(rangeMask(1:end-1) ~= rangeMask(2:end));



%% EXTRACT TIME RANGE

% Note the time input doesn't need to correspond to an entry in adcpInput.t. So
% we need to allow for this and can't just use find(). Times in adcpInput.t must
% be strictly increasing, however, which is useful.
tInds = 1:numel(adcpInput.t);

% Some horrid, ugly code to find the nearest index in adcpInput to the input
% time.
i1 = find(adcpInput.t <= t, 1, 'first');
i2 = find(adcpInput.t >  t, 1, 'last' );
if isempty(i1)
    timeInd = i2;
elseif isempty(i2)
    timeInd = i1;
elseif adcpInput.t(i1) == t
    timeInd = i1;
else % find the nearest
    t1 = adcpInput.t(i1);
    t2 = adcpInput.t(i2);
    [~, whichOne] = min([(t-t1), (t2-t)]);
    if whichOne == 1
        timeInd = i1;
    else
        timeInd = i2;
    end
end

% Now get the switch points either side of the time of interest
lowerSwitch = find(switchInds <= timeInd, 1, 'last' );
upperSwitch = find(switchInds >  timeInd, 1, 'first');

% And get the window out of the ADCP data structure
[adcpWindow] = getWindow(adcpInput, 'indexRange', [lowerSwitch upperSwitch]);



%% GRAPHICAL OUTPUT

% Get the index range to plot (50% of the half cycle again either way. Use a
% monotonic assumption to do this, since it's easier and doesn't affect results
lowerPlot = ceil(lowerSwitch - (upperSwitch-lowerSwitch)*0.5);
upperPlot = floor(upperSwitch + (upperSwitch-lowerSwitch)*0.5);
if lowerPlot <= 0
    lowerPlot = 1;
end
if upperPlot > numel(adcpInput.t)
    upperPlot = numel(adcpInput.t);
end
plotMask = false(size(adcpInput.t));
plotMask(lowerPlot:upperPlot) = true;



% Render the figure
raiseFigure('Half Cycle Velocities')
clf
subplot(2,2,1)
plot(adcpInput.t(plotMask), u(plotMask), 'b-')
xlabel('u (m/s)')
ylabel('Time [datenum form] (s)')
hold on
plot([adcpInput.t(lowerSwitch) adcpInput.t(upperSwitch)], [0 0], 'ro')

subplot(2,2,2)
plot(adcpInput.t(plotMask), v(plotMask), 'g-')
xlabel('v (m/s)')
ylabel('Time [datenum form] (s)')
hold on
plot([adcpInput.t(lowerSwitch) adcpInput.t(upperSwitch)], [0 0], 'ro')

subplot(2,2,3)
plot(adcpInput.t(plotMask), w(plotMask), 'm-')
xlabel('w (m/s)')
ylabel('Time [datenum form] (s)')
hold on

subplot(2,2,4)
plot(adcpInput.t(plotMask), degreesNautical(plotMask), 'k-')
xlabel('Flow Direction  (Nautical \circ)')
ylabel('Time [datenum form] (s)')
hold on

if ~strcmp(opts.filter,'none')
    
    subplot(2,2,1)
    plot(adcpInput.t(plotMask), uFilt(plotMask), 'r--')
    hold on

    subplot(2,2,2)
    plot(adcpInput.t(plotMask), vFilt(plotMask), 'r--')
    hold on
    
end


end % end main function