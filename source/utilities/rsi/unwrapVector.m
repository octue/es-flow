function [out] = unwrapVector(in, ambiguityVelocity, varargin)
%UNWRAPVECTOR Unwraps a nortek vector based on the criterion that a jump of more
%than the ambiguity velocity in between two consecutive samples consititues a
%wrap around
% Syntax:
%       [uOut] = unwrapVector(uIn, ambigVelocity) accepts beam coordinate
%       inputs and returns beam coordinate outputs unwrapped by the ambiguity
%       velocity
%
%       [uOut] = FunctionName(uIn, ambigVelocity, T) accepts cartesian
%       coordinate velocities as input and returns cartesian coordinate outputs
%       unwrapped by the ambiguity velocity in the beam coordiante system
%
%
% Inputs:
%
%       uIn             [n x 3]     Input velocity components in m/s in beam
%                                   coordinates (or in cartesian coords if T
%                                   also passed)
%
%      ambigVelocity    [1 x 1]     Ambiguity velocity in m/s. Allowable values
%                                   and ping ranges to correspond:
%
%                 Nominal Velocity Range (m/s)    Ping Interval (sec)
%                 0.03                            5.0e-4
%                 0.10                            3.0e-4
%                 0.30                            1.4e-4
%                 1.00                            7.0e-5
%                 2.50                            4.0e-5
%                 4.00                            2.5e-5
%
%       T               [3 x 3]     Transformation matrix from beam to cartesian
%                                   coords. Given in the vectrino file header,
%                                   e.g.
%                                     T*[b1; b2; b3] = [u; v; w]
%                                     inv(T)*[u; v; w] = [b1; b2; b3]
%   
% Outputs:
%
%       uOut            [n x 3]     Output velocity components in m/s in beam
%                                   coordinates (or in cartesian coords if T
%                                   also passed)
%
% Future improvements
%
%   [1] Use of a wavelet filter to provide improved unwrapping performance in
%       the presence of significant fluctuations in the data.
%
% References:
%
%   [1] Nortek website discussion page http://www.nortek-as.com/en/knowledge-center/forum/general-comments/530960035
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
% Revision History:        	23 February 2015    Created
%
% Copyright (c) 2015 Ocean Array Systems, All Rights Reserved.
%
% Sort variable arguments to give transformation matrices
if nargin > 2
    T = varargin{1};
    invT = inv(T);
else
    T = eye(3);
    invT = eye(3);
end
    
% Throw error if ambiguity velocity is nonstandard for the instrument
if ~any(bsxfun(@eq, ambiguityVelocity, [0.03 0.1 0.3 1 2.5 4]))
    error('Ambiguity velocity value not one of the predefined set.')
end

% Transform in if necessary
in = invT*in';

% Scale everything so the ambiguity velocity is = pi, and the wraparound is 2pi
in = in*pi/ambiguityVelocity;

error('this doesnt work - check againts real data')
% Unwrap
out(1,:) = unwrap(in(1,:),pi/2);
out(2,:) = unwrap(in(2,:),pi/2);
out(3,:) = unwrap(in(3,:),pi/2);

% Transform back
out = ambiguityVelocity*(T*out)'/pi;
end

