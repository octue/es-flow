function [tf] = isMonotonic(adcpData, varargin)
%ISMONOTONIC Returns true if input adcp data series is monotonic in time 
%
% Syntax:  
%       [tf] = isMonotonic(adcpInput) Returns a boolean value which is true if
%       the timebase of the input adcpData structure is monotonic (uniform time
%       increment between each profile, within a margin of error corresponding
%       to 10*eps('single')
%
%       [tf] = isMonotonic(adcpInput, epsValue) Allows user to specify a
%       threshold other than the default 10*eps('single'), below which
%       differences in successive timebase values are considered to be
%       negligible
%
%       [tf] = isMonotonic(adcpInput, epsValue, 'seconds') Allows user to specif
%       a threshold in units of seconds, rather than in the MATLAB datenum base
%       
% Inputs:
%
%       adcpInput       structure       An OAS standard adcp data structure.
%                                       See help('loadADCP') for fields
%                                       reference.
%
%       epsValue        [1 x 1]         The value below which variation between
%                                       successive entries in the timebase
%                                       vector is considered noise. Default
%                                       10*eps('single') (equivalent to approx
%                                       0.103s error for values stored in
%                                       datenum form).
%
% Outputs:
%
%   	tf              [1 x 1]         Returns true if the timebase is
%                                       monotonic within the error constraints
%                                       set, returns false otherwise.
%
% Future Improvements:      
%
%       [1] Include ability to force monotonic timebase by rounding to a
%           specific decimal place. 
%
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
% Created:          25 July 2014
% Revisions:        

% Checks on input size, class and type

% Set up for different values
if nargin > 2
    epsValue = datenum([0 0 0 0 0 varargin{1}]);
    
elseif nargin > 1
    % Use input epsValue in dateNum format
    epsValue = varargin{1};
    
else
    % Use default epsValue
    epsValue = double(10*eps('single'));
end

% Check for monotonic times
tf = true;
dt = diff(adcpData.t);
if any((dt(2:end) - dt(1:end-1)) > epsValue)
    warning('MATLAB:NonMonotonic','isMonotonic: Non monotonic data found in input file series.')
    tf = false;
end


