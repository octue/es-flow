function [tauStarA, tauStarB] = getShear(eta, Pi, S, zeta, beta)
%GETSHEAR gets normalised shear stress profiles as a function of eta
%
% Syntax:  
%       [tauStarA, tauStarB] = getShear(eta, Pi)
%
% Inputs:
%
%       eta         [nZ x 1]    Nondimensional height z/deltac above the wall at
%                               which Reynolds Stresses are required. Eta values
%                               must ascend but not necessarily be monotonic.
%       
%   	Pi          [1 x 1]     Coles wake parameter Pi
%
% Outputs:
%
%   	tauStarA    [nZ x 1]    Shear stress profile using Lewkowicz (1982)
%                               formulation. This is the sink flow contribution
%                               of the Type A eddies alone.
%
%       tauStarB    [nZ x 1]    Contribution of Type B eddies to the overall
%                               shear stress.
%
% References:
%
%   [1] Perry AE and Marusic I (1995) A wall-wake model for turbulent boundary
%       layers. Part 1. Extension of the attached eddy hypothesis J Fluid Mech
%       vol 298 pp 361-388
%
% Future Improvements:  
%
%   [1] Optional different mean profile formulations including account for the
%       free surface 
%
%   [2] Support directional variation with height; i.e. compatible with mean
%       profile formulations using U(y) 
%
%   [3] Added formulation for contribution of smaller Type C eddies
%
% Author:                   V. Gupta and T. Clark
% Work address:             Ocean Array Systems Ltd
%                           Hauser Forum
%                           3 Charles Babbage Road
%                           Cambridge
%                           CB3 0GT
% Email:                    vikrant.gupta@oceanarraysystems.com
% Website:                  www.oceanarraysystems.com
%
% Created:                  31 July 2014        Created
%                           19 April 2015       Modified by T. Clark from the
%                                               original function
%                                               u1u3_fromParams() into getShear.
%                                               Slightly modified the quilibrium
%                                               sink flow formulation for
%                                               consistency with the Lewkowicz
%                                               model (see Perry and Marusic eq.
%                                               51). Altered outputs to be
%                                               nondimensional.
%
% Copyright (c) 2014-2015 Ocean Array Systems, All Rights Reserved.



% Get the gradient of the velocity deficit wrt lambda. See hand derivation PDF
% attached to source code (nb we multiply by dEta/dLambdaE so derivative of the
% velocity deficit is with respect to LambdaE not wrt Eta
dUds_dEta = -1/(kappa*eta) + (Pi/kappa)*( eta*(12-(2/Pi)) + (eta.^2)*((6/Pi)-12) + (eta.^3)*(-5/Pi));
dUds_dLambdaE = -1*exp(-1*lambdaE).*fliplr(dUds_dEta);



