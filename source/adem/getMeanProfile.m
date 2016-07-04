function [Ux] = getMeanProfile(Pi, S, deltac, U1, z)
%GETMEANPROFILE Determine mean boundary layer velocity u(z) 
% Uses the using the Lewkowicz (1982) formulation (Perry and Marusic eq.9) of
% the Coles wake function to determine U(z) from parameters Pi, S, deltac, U1.
%
% Syntax:  
%       [Ux] = getMeanProfile(Pi, S, deltac, U1, z)
%       Determines u(z) using the approach described in Ref. 1.
%
% Inputs:
%       z           [nZ x 1]    Height in m above the wall at which boundary
%                               layer profile is required.
%       
%   	Pi          [1 x 1]     Coles wake parameter Pi
%
%       S           [1 x 1]     Ratio between free stream and friction velocity
%                               S = U1/Utau
%
%       deltac      [1 x 1]     The boundary layer thickness in m
%
%       U1          [1 x 1]     The free stream speed in m/s
%
% Outputs:
%
%       Ux          [nZ x 1]    Streamwise velocity in m/s at points
%                               corresponding to the heights in z
%
%       dUds_dEta    [nZ x 1]   Gradient of the normalised velocity deficit
%                               Uds = (U1-Ux)/Utau with respect to 
%                               the eta coordinate.
%
% References:
%
%   [1] Perry AE and Marusic I (1995) A wall-wake model for turbulent boundary
%       layers. Part 1. Extension of the attached eddy hypothesis J Fluid Mech
%       vol 298 pp 361-388
%
% Future Improvements:  
%
%   [1] Constrained but variable value of kappa (presently kappa = 0.41 always)
%
%   [2] Ability to use different parameter models (e.g. use a
%       structure (or an object?) to give input parameters
%
%   [3] Support directional variation with height; i.e. accept U(y)
%
% Other m-files required:   none
% Subfunctions:             none
% Nested functions:         none
% MAT-files required:       none
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
% Revisions:        07 April 2015       Created
%                   18 April 2015       Properly documented with OAS header,
%                                       minor change of output name from u to
%                                       Ux. Renamed getMean to getMeanProfile.
%
% Copyright (c) 2014-2015 Ocean Array Systems, All Rights Reserved.

% Define von Karman constant
kappa = 0.41;

% S is the ratio between free stream velocity and the wall friction velocity -
% use it to get U1 (eq.7 Perry and Marusic)
Utau = U1/S;

% Nondimensional wall distance
eta = z/deltac;

% Coles wake function, as a function of Pi and eta, using the Lewkowicz (1982)
% formulation (Perry and Marusic eq.9)
Wc = 2*eta.^2.*(3-2*eta) - eta.^2.*(1-eta).*(1-2*eta)/Pi;

% Determine the analytical profile of the boundary layer
uDeficit = -log(eta)/kappa + (Pi/kappa)*2 - (Pi/kappa)*Wc;
u = U1 - uDeficit*Utau;
uDeficit1 = -log(1)/kappa + (Pi/kappa)*2 - (Pi/kappa)*2;

% Determine the fixed value outside the boundary layer
uBar1 = U1 - uDeficit1*Utau;
u(eta>=1) = uBar1;

% Outputs
Ux = u;


end
