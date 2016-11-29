function [windows] = windowMR(tYD, uvw, heading, pressure, temperature, height, varargin)
%WINDOWVECTOR Takes vector data and split it into windows
%
% Syntax:
%       [windows] = windowVector(tYD, uvw, heading, height, waterDepth) Extracts
%       segments from the input vector data consisting of 10 minute long windows
%       of flow.
%
%       [windows] = windowVector(..., 'TimeOffset', timeOffset) 
%       Extracts segments from the input vector data consisting of 10 minute
%       long windows of flow whose start points are offset by offsetTime
%       seconds. Note that the nearest offset time creating noninteger indices
%       is selected (exact offset can be determined by comparing output time
%       fields)
%
%
% Inputs:
%
%       tYD             [n x 1]     Time in year day form of each entry
%
%       uvw             [n x 3]     Velocity U,V,W wrt the nemo reference frame
%                                   at time tYD
%
%       heading         [n x 1]     Heading of the NEMO in degrees at time tYD
%
%       temperature     [n x 1]     Temperature at time tYD
%
%       height          [n x 1]     Height in m from the seabed at time tYD
%
%       pressure        [n x 1]     Pressure in Pascals at time tYD
%
%       timeOffset      [n x 1]     Offset in seconds between starts of
%                                   successive bins. This allows bins to
%                                   overlap. If left empty bins are consecutive
%
% Outputs:
%
%       windows          [m x 1] structure
%                                   Structure for m windows containing the
%                                   following fields:
%
%        .startTimeYD   [1 x 1]     Start time in year day format.
%
%        .uvw           [p x 3]     Velocity extracted either directly from uvw
%                                   input or from a filtered version of the
%                                   input (depending on variable arguments) for
%                                   a 10 minute window.
%
%        .uvwMean       [1 x 3]     Mean uvw in 10 minute window
%
%        .time          [p x 1]     Time in seconds from startTimeYD
%
%        .waterDepth    [1 x 1]     Mean water depth during the window
%
%        .height        [1 x 1]     Mean NEMO height from the seabed in m during
%                                   the window.
%
%        .heading       [1 x 1]     Mean heading in degrees during the
%                                   window.
%
% See Also: WP1
%
% References:
%
%   [1] T. H. Clark "Measurement of three-dimensional coherent fluid structure
%       in high Reynolds Number turbulent boundary layers. University of
%       Cambridge, February 2012
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
% Revision History:        	23 February 2015
%
% Copyright (c) 2015 Ocean Array Systems, All Rights Reserved.

% Set nondefault parsed options
opts.TimeOffset = 10*60; % 10 minutes
opts = parse_pv_pairs(opts,varargin);

% Initialise the structure for output
windows = struct('startTimeYD', [],  'time', [], 'uvw', [], 'heading', [],'temperature',[], 'pressure', [], 'height', []);

% Get bin start indices
timeOffsetYD = opts.TimeOffset/(60*60*24);
startTimes = tYD(1):timeOffsetYD:tYD(end);

% Attempt to add all the bins into a structure. NB if the windows aren't
% contiguous, or are out of bounds, valid is false and the bin isn't added. Thus
% we don't know a priori what size the output is.
validCtr = 1;
for i = 1:numel(startTimes)
    [valid, window] = getWindow(startTimes(i), tYD, uvw, heading, temperature, pressure, height);
    if valid
        windows(validCtr) = window;
        validCtr = validCtr + 1;
    else
        disp(['Window ' num2str(i) ' invalid due to discontinuous timeseries'])
    end
end

end % end main function

function [valid, window] = getWindow(startTimeYD, tYD, uvwFilt, heading, temperature, pressure, height)
    
    
% Get endTime, 10 minutes (600 seconds) after the startTime, converted to decimal days
endTimeYD = startTimeYD + (10*60)/(60*60*24);

% Get start and end indices of the window
startInd = find(tYD >= startTimeYD, 1, 'first');
endInd = find(tYD <= endTimeYD, 1, 'last');

% Return an invalid entry if this puts us out of bounds (i.e. empty matrices are
% returned for one or both)
if isempty(startInd) || isempty(endInd)
    valid = false;
    window = [];
    return
end

try

    % Get the window out and assign to the structure
    window.startTimeYD = tYD(startInd);
    window.time = (tYD(startInd:endInd) - window.startTimeYD)*60*60*24;
    window.uvw = uvwFilt(startInd:endInd,:);
    window.heading = heading(startInd:endInd);
    window.temperature = temperature(startInd:endInd);
    window.pressure = pressure(startInd:endInd);
    window.height = height(startInd:endInd);
    valid = true;
    
    % Throw an error if the timeseries varies by more than 1% per timestep (this
    % indicates discontinuities in data, e.g. instrument being switched off then
    % on)
    dt = diff(window.time);
    dtVariation = dt./dt(1);
    if any(dtVariation > 1.01) || any(dtVariation < 0.99)
        error('windowVector:DiscontinuousTime','Nonuniform timeseries in window')
    end
    
catch
    
    % If a dicontinuity has been detected, catch the error and return an invalid
    % window
    valid = false;
    window = [];
    
end

end % end subfunction getWindow()





