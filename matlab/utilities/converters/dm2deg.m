function [deg, varargout] = dm2deg(d, m, varargin)
%DM2DEG Converts a degrees-minutes pair or pais to a decimal degrees expression
%
% Syntax:
%       [deg] = FunctionName(d,m) Converts degrees minutes to decimal degrees
%
%       [deg, deg2] = FunctionName(d,m,d2,m2) Allows an additional set of pairs
%       to be converted (e.g. so latitude and longitude can be converted in one
%       line)
%
%
% Inputs:
%
%       d          [n x 1]     Degrees
%
%       m          [n x 1]     Minutes in the range 0= < m < 60 
%
% Outputs:
%
%       deg         [n x 1]     Decimal degrees expression of [d m].
%
%       deg2         [n x 1]     Decimal degrees expression of [d2 m2].
%
% See Also: DMS2DEG, deg2utm, utm2deg, deg2dm, deg2dms
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
% Revision History:        	02 April 2015
%
% Copyright (c) 2014-2015 Ocean Array Systems, All Rights Reserved.


deg = d+m/60;

if ((nargout < 2) && (nargin > 2)) || ((nargout < 2) && (nargin > 2))
    error('Invalid number of inputs or outputs')
end

if nargin > 2
    
    varargout{1} = varargin{1} + varargin{2}/60;
    
end
    
    
end

