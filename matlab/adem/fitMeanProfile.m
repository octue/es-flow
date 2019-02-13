function profile = fitMeanProfile(z, Ux, weights, x0, constrained, type)
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
%       constrained [4 x 1]     logical mask. Where true, the corresponding
%                               value of x0 is constrained to be a solution. 
%                               Where false, the value is allowed to vary.
%
%       type        char        string denoting the model type used. Currently
%                               only Lewkowicz 1982 model implemented.
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
% Work address:             Octue
%                           Hauser Forum
%                           3 Charles Babbage Road
%                           Cambridge
%                           CB3 0GT
% Email:                    tom@octue.com
% Website:                  www.octue.com
%
% Copyright (c) 2014-2018 Octue, All Rights Reserved.

% TODO Input checks on data type and dimension

% Adding a switch based on type which is currently redundant
if nargin < 6
    type = 'lewkowicz';
end

% The input z is a 1 x 40 double precision array, while Ux is a [39 x 1] double
% precision array

% Force inputs to column vectors and check size
Ux = Ux(:);
z = z(:);
if any(size(Ux) ~= size(z))
    error('Inputs Ux and z must be of the same size')
end

% Define the boundary layer data that we know (mean profiles).
xData = z(~isnan(Ux));
yData = Ux(~isnan(Ux));

% Ensure double precision inputs
% TODO Ensure that sensible precision data is saved in the .adcp files and
% manage correct precision throughout.
xData = double(xData);
x0 = double(x0);
yData = double(yData);

% Handle an entirely NaN input
if isempty(yData)
    profile.Pi      = NaN;
    profile.S       = NaN;
    profile.deltac  = NaN;
    profile.U1      = NaN;
    profile.Utau    = NaN;
    profile.kappa   = NaN;
    profile.resnorm = NaN;
    profile.type    = type;
    return
end

% Set weightings and mask out nan entres
if isempty(weights)
    weights = ones(size(z));
else
    warning('weighting not implemented')
end
weights = weights(~isnan(Ux));

% Set constraints, ensure logical input, remove constrained values from x0
if nargin < 5 || isempty(constrained)
    constrained = false(size(x0));
else
    % Ensure input is logical
    constrained = logical(constrained);
end
cons = nan(size(constrained));
cons(constrained) = x0(constrained);
x0 = x0(~constrained);

% Define von Karman constant
kappa = 0.41;

% Suppress display of optimisation progress
opts = optimoptions('lsqcurvefit','Display','iter','Algorithm', 'levenberg-marquardt')

% Run the curve fit, for different types of boundary layer profile. Using
% lsqcurvefit is a more robust approach than the previous minimisation approach,
% which used a Nelder-Meade simplex search and was quite flaky.
switch type
    
    case 'lewkowicz'
        try
        [xFit, profile.resnorm] = lsqcurvefit(@(x,xData) lewkowicz(x, xData, kappa, weights, cons), x0, xData, yData, [], [], opts);
        catch me
            save('debug.mat','xData','yData','opts','x0','kappa','weights')
            throw(me)
        end
        
        % Splice fitted values in xFit with constrained values to make a
        % parameters vector, then output results in a structure
        params = cons;
        params(isnan(cons)) = xFit;
        profile.Pi      = params(1);
        profile.S       = params(2);
        profile.deltac  = params(3);
        profile.U1      = params(4);
        profile.Utau    = profile.U1/profile.S;
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

% Handle any complex output
if ~isreal(profile.Pi)
    save('debug_complex.mat','xData','yData','opts','x0','kappa','weights')
    profile.Pi      = NaN;
    profile.S       = NaN;
    profile.deltac  = NaN;
    profile.U1      = NaN;
    profile.Utau    = NaN;
    profile.kappa   = NaN;
    profile.resnorm = NaN;
    return
end

end



function [fx] = lewkowicz(x, xData, kappa, weights, cons)

% Splice varying values in x with constrained values to make a parameters vector
params = cons;
params(isnan(cons)) = x;

% Extract as named variables for code clarity
Pi      = params(1);
S       = params(2);
deltac  = params(3);
U1      = params(4);

% Get parameters out of xData
z = xData(:,1);

% Nondimensional wall distance
eta = z/deltac;
% S is the ratio between free stream velocity and the wall friction velocity -
% use it to get U1 (eq.7 Perry and Marusic)
Utau = U1/S;

% sprintf('pi %f', Pi)
% sprintf('u_inf %f', U1)
% sprintf('u_tau %f', Utau)
% sprintf('delta_c %f', deltac)

% Coles wake function, as a function of Pi and eta, using the Lewkowicz (1982)
% formulation (Perry and Marusic eq.9)
Wc = 2*eta.^2.*(3-2*eta) - eta.^2.*(1-eta).*(1-2*eta)/Pi;

% Determine the analytical profile of the boundary layer
uDeficit = -log(eta)/kappa + (Pi/kappa)*2 - (Pi/kappa)*Wc;
u = -uDeficit*Utau + U1;
uDeficit1 = -log(1)/kappa + (Pi/kappa)*2 - (Pi/kappa)*2;
uBar1 = -uDeficit1*Utau + U1;

% Treat values outside the boundary layer as constant velocity
u(eta>=1) = uBar1;

% Weight according to the confidence in uBar throughout the water column.
% TODO introduce some kind of weighting - this doesn't work!
fx = u;

end
