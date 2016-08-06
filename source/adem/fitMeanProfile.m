function profile = fitMeanProfile(z, Ux, weights, x0, type)
%FITMEANPROFILE Least square fit to find mean boundary layer parameters.
% Perform a least squares fit to ascertain mean boundary layer parameters, for a
% range of different analytical profiles.
%
% Profiles available:
%   'lewkowicz'
%       Returns Pi, S, deltac, U1, kappa, U_tau for the experimentally
%       determined 2d boundary layer f(Ux,z). Uses the using the Lewkowicz
%       (1982) formulation (Perry and Marusic eq.9) of the Coles wake function.
%
%   'logarithmic'
%       Not implemented yet.
%
%   'exponential'
%       Not implemented yet.
%
% Syntax:  
%       profile = fitMeanProfile(z,Ux, weight, x0)
%       Determines parameters describing an equilibrium or quasi-equilibriium
%       boundary layer for use with analyses in Ref [1].
%
% Inputs:
%       z           [nZ x 1]    Height in m above the wall at which boundary
%                               layer profile was measured. Ascending but not
%                               strictly monotonic.
%       
%       Ux          [nZ x 1]    Streamwise velocity in m/s at points
%                               corresponding to the heights in z
%
%       weights     [nZ x 1]    Relative confidence weightings ranging between 0
%                               and 1 for the Ux measurements at corresponding
%                               heights. 0 = no confidence (likely just noise)
%                               while 1 = best confidence. Leave empty for no
%                               weightings (set to all 1s).
%
%       x0          [4 x 1]     Initial guesses for [Pi, S0, deltac0, U10]. 
%
% Outputs for type lewkowicz:
%
%       profile     struct with the following fields:
% 
%           .Pi          [1 x 1]    Coles wake parameter Pi
% 
%           .S           [1 x 1]    Ratio between free stream and friction
%                                   velocity 
%                                   S = U1/Utau
% 
%           .deltac      [1 x 1]    The boundary layer thickness in m
% 
%           .U1          [1 x 1]    The free stream speed in m/s
% 
%           .Utau        [1 x 1]    The skin friction velocity in m/s
% 
%           .kappa       [1 x 1]    Kappa, the von karman constant used for the
%                                   fit. Presently 0.41 always (see future
%                                   improvements)
%
%           .resnorm     [1 x 1]    Normalised residual associated with the fit.
%
%           .type        string     Type of fit, 'lewkowicz'.
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
%   [2] Input value of deltac (reduces dimensionality of the fit)
%
%   [3] Automated first guess based on calculation from the profile input
%
%   [4] Optional different fits including account for the free surface
%
%   [5] Revised output to account for the different parameter models (e.g. use a
%       structure (or an object?) to contain result
%
%   [6] Add the fit quality to the output so we have an additional metric which
%       can be used in smoothing
%
%   [7] Support directional variation with height; i.e. accept U(y)
%
%   [8] Add option to weight the confidence in the input vector
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
% Copyright (c) 2014-2015 Ocean Array Systems, All Rights Reserved.

% TODO Input checks on data type and dimension

% Adding a switch based on type which is currently redundant
if nargin < 5
    type = 'lewkowicz';
end

% Define the boundary layer data that we know (mean profiles).
xData = z;
yData = Ux;

% Set weightings if empty is given
if isempty(weights)
    weights = ones(size(z));
end

% Define von Karman constant
kappa = 0.41;

% Suppress display of optimisation progress
opts = optimoptions('lsqcurvefit','Display','off');

% Run the curve fit, for different types of boundary layer profile. Using
% lsqcurvefit is a more robust approach than the previous minimisation approach,
% which used a Nelder-Meade simplex search and was quite flaky.
switch type
    
    case 'lewkowicz'
        [xFit, profile.resnorm] = lsqcurvefit(@(x,xData) lewkowicz(x, xData, kappa, weights), x0, xData, yData, [],[],opts);
        profile.Pi      = xFit(1);
        profile.S       = xFit(2);
        profile.deltac  = xFit(3);
        profile.U1      = xFit(4);
        profile.Utau    = U1/S;
        profile.kappa   = kappa;
        
    case 'logarithmic'
        error('not implemented yet')
        
    case 'exponential'
        error('not implemented yet')
        
    otherwise
        error('Unknown velocity profile type. Try ''lewkowicz'', ''logarithmic'', or ''exponential''.')
        
end

% Outputs
profile.type = type;

end



function [fx] = lewkowicz(x, xData, kappa, weights)

% Extract for clarity
Pi      = x(1);
S       = x(2);
deltac  = x(3);
U1      = x(4);

% Get parameters out of xData
z = xData(:,1);

% Nondimensional wall distance
eta = z/deltac;

% S is the ratio between free stream velocity and the wall friction velocity -
% use it to get U1 (eq.7 Perry and Marusic)
Utau = U1/S;

% Coles wake function, as a function of Pi and eta, using the Lewkowicz (1982)
% formulation (Perry and Marusic eq.9)
Wc = 2*eta.^2.*(3-2*eta) - eta.^2.*(1-eta).*(1-2*eta)/Pi;

% Determine the analytical profile of the boundary layer
uDeficit = -log(eta)/kappa + (Pi/kappa)*2 - (Pi/kappa)*Wc;
u = -uDeficit*Utau + U1;
uDeficit1 = -log(1)/kappa + (Pi/kappa)*2 - (Pi/kappa)*2;
uBar1 = -uDeficit1*Utau + U1;
u(eta>=1) = uBar1;

% Weight according to the confidence in uBar throughout the water column.
fx = weights.*u;

end
