function dispnow(varargin)
%DISPNOW Displays output to the command line with a datestamp, as per the
%requirement for TurbineGrid outputs. Now accepts sprintf style multiple
%arguments, so use it like sprintf.
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

disp([datestr(now) '    ' sprintf(varargin{:})])