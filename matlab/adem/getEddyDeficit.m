function [h, lambda] = getEddyDeficit(type, varargin)
%GETEDDYDEFICIT Returns eddy deficit signature function h
% This represents the contribution to the gradient of the velocity deficit that
% an individual structure will contribute to a boundary layer.
%
% Syntax:
%       [h, lambda] = getEddyDeficit('type')
%       Returns eddy deficit function h for the eddy type described in
%       the input string, presently 'A', 'B1', 'B2', 'B3', 'B4' which are the
%       eddies described in Ref 1.
%
%       [h, lambda] = getEddyDeficit('type', lambda)
%       Returns eddy deficit function h for prescribed wall coordinates lambda
% Inputs:
%
%       type            string      The eddy type for which to calculate J.
%                                   Current accepted values are  'A', 'B1',
%                                   'B2', 'B3', 'B4', which correspond with the
%                                   eddy types in Ref [1].
%
%       lambda          [nZ x 1]    The remapped wall coordinate at which J is
%                                   given.
%
% Outputs:
%
%       h               [nZ x 1]    The contribution to the mean shear made by
%                                   an eddy at each wall normal point lambda
%
%       lambda          [nZ x 1]    The remapped wall coordinate at which J is
%                                   given.
%
% See Also: GETEDDYINTENSITY.M
%
% Future Improvements:
%
%   [1] Computation should be done analytically (for a given vortex core
%       function like rankine) rather than assuming line vortex.
%
%   [2] Possible additional modification to take into account a free surface
%       image or additional structure types
%
% References:
%
%   [1] Perry AE and Marusic I (1995) A wall-wake model for turbulent boundary
%       layers. Part 1. Extension of the attached eddy hypothesis J Fluid Mech
%       vol 298 pp 361-388
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
% Revision History:        	19 April 2015       Created from getEddyIntensity.m
%
% Copyright (c) 2014-2015 Ocean Array Systems, All Rights Reserved.


%% VERTICAL POSITIONS

% And infer J which is equal; just remapped to a different space
% Compute the z/delta spacing from lambda the logarithmic spacing
if nargin > 1
    lambda = varargin{1};
else
    lambda = linspace(0, log(1/0.001), 1000);
end
eta = 1./exp(lambda);

% Initialise f
f = zeros(size(eta));

% Switch according to type requested; and determine f(z/delta) using figure 3 in
% Perry and Marusic
switch lower(type(1))
    case 'a'
        
        mask = (eta>=0) & (eta <= 1);
        f(mask) = 1;
        
    case 'b'
        
        mask = (eta>=0.5) & (eta <= 1);
        f(mask) = 1;
                    
    otherwise
        
        error('MATLAB:InvalidInput','Unrecognised eddy type string input. See help(''getEddyIntensity'') for valid eddy types')
        
end

% h(lambda) = f(z/delta):
h = f;

