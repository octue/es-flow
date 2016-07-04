function [varargout] = timeYearDay(varargin)
%TIMEYEARDAY Convertstimes between datevec and year-day format
%
% Syntax:
%       [decDays, year] = timeYearDay(dv) converts datevec to decimal days form
%
%       [dv] = timeYearDay(ddays, year) converts decimal days to datevec form
%
%
% Inputs:
%
%       dv              [m x 6]     datevecs (see help('datevec'))
%
% Outputs:
%
%       decDays         [m x 1]     Decimal days elapsed since start of year
%
%       year            [m x 1]     Year (equal to datevec(:,1) )
%
% See Also: DATEVEC.M
%
% Future Improvements:
%
%   [1] Full checks. Currently implemented to be compatible with RSI's decDays
%       form, but I don't think theirs is right - check what happens e.g. in the
%       early hours of january 1st.
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
% Revision History:        	04 March 2015       Created
%
% Copyright (c) 2015 Ocean Array Systems, All Rights Reserved.

% Cumulative days in each month
monthDays = cumsum([0 31 28 31 30 31 30 31 31 30 31 30 31]);

switch nargin
    case 1
        
        % CONVERT DATEVEC TO DECIMAL DAYS
        dv = varargin{1};
        year = dv(:,1);                 % Year
        ddays = monthDays(dv(:,2))';            % Months
        ddays = ddays + (dv(:,3));              % Days - I think that this should be dv(:,3)-1 but sticking with the RSI convention for now
        ddays = ddays + (dv(:,4)/24);           % Hours
        ddays = ddays + (dv(:,5)/(24*60));      % Minutes
        ddays = ddays + (dv(:,6)/(24*60*60));   % Seconds
        
        varargout{1} = ddays;
        varargout{2} = year;
        
    case 2
        
        % CONVERT DECIMAL DAYS TO DATEVEC
        
        ddays = varargin{1};
        years = varargin{2};
        
        dv = [years, zeros(numel(ddays),5)];
        
        for i = 1:numel(years)
            t = datetime([years(i) 1 0 0 0 0]);
            t = t + ddays(i);
            dv(i,:) = datevec(t);
        end
        
        varargout{1} = dv;

end

