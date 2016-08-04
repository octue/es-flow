function [adcpData] = correctTime(adcpData, drift, varargin)
%CORRECTTIME Corrects for drift and forces strict monotonicity of ADCP timebase. 
%
% Syntax:  
%       [tf] = correctTime(adcpInput, drift) Returns an ADCP data structure
%       whose timebase is corrected for an externally measured drift over a
%       deployment period. Once corrected, the timebase is forced to a strictly
%       monotonic series. Inputting timeseries whose time increments are not
%       approximately monotonic (i.e. successive differences in time increments
%       vary by less than a low pre-set threshold) produces an error to protect
%       against re-timebasing data not suited for this function.
%
%       [tf] = correctTime(adcpInput, drift, epsValue) Allows user to specify a
%       threshold other than the default 10*eps('single'), below which
%       differences in successive timebase increments are considered to be
%       negligible.
%
%       [tf] = correctTime(adcpInput, drift, epsValue, 'seconds') Allows user to
%       specify a non-default threshold in units of seconds, rather than in the
%       MATLAB datenum form.
%       
% Inputs:
%
%       adcpInput       structure       An OAS standard adcp data structure.
%                                       See help('loadADCP') for fields
%                                       reference.
%
%       drift           [1 x 1]         Drift in seconds relative to an external
%                                       timebase. Positive drift suggests the
%                                       ADCP was clocking faster than the
%                                       external source and therefore had an
%                                       apparently greater time elapsed during
%                                       a deployment.
%
%       epsValue        [1 x 1]         The value, in percent, above which
%                                       variation in sampling period between
%                                       successive samples is considered to
%                                       represent a discontinuity in the data,
%                                       rather than normal variation due to
%                                       clocking rates, drift, temperature
%                                       fluctuations and roundoff errors.
%                                       Default is 1%. Set to Inf to ignore
%                                       variations (not advised).
%
% Outputs:
%
%   	adcpData        structure       An OAS standard ADCP data structure,
%                                       with the timebase refactored to account
%                                       for drift and forced to a strictly
%                                       monotonic timeseries.
%
% Future Improvements: 
%
%       [1] Ability to resample data to an externally imposed timebase.
%
%       [2] Addition of a starting drift to account for drift occuring between
%           the synchronisation with external source and the first ping (assumed
%           negligible at present).
%
%       [3] Incorporation of drift rates to make the usage a bit more versatile.
%
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
% Revision History:        	25 July 2014        Created
%                           06 April 2015       Altered header to conform to OAS
%                                               standard and added a future
%                                               improvement request.
%
% Copyright (c) 2014-2015 Ocean Array Systems, All Rights Reserved.

% Checks on input size, class and type

% Set up for different values
if nargin > 2
    % Use input epsValue
    epsValue = varargin{1}/100;
else
    % Use default epsValue in seconds
    epsValue = 0.01;
end


% Get the approximate timebase
startTime = adcpData.t(1);
approxEndTime = adcpData.t(end);

% Subtract the drift (which is given in seconds, positive when clocking faster
% than the external sync source) to correct the timebase
correctedEndTime = addtodate(approxEndTime, -1*drift, 'second');

% Return an error if the timebase varies by more than 1% step to step.
dt = diff(adcpData.t);
if any(((dt(2:end) - dt(1:end-1))./dt(1:end-1)) >= epsValue)
    error('MATLAB:InvalidTimebase', ['Sampling period dt fluctuates by more than ' num2str(epsValue*100) ' percent between adjacent samples'])
end

% Impose strictly monotonic timebase between the start and the end times
adcpData.t = linspace(startTime, correctedEndTime, numel(adcpData.t));

    
