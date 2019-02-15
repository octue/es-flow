function [r13A, r13B] = getReynoldsStress13(eta, Pi, S, zeta, beta)
%GETREYNOLDSSTRESS13 Gets Reynolds Stress profiles due to type A and B eddies
% Adopts the approach of Perry and Marusic 1995 and uses the using the Lewkowicz
% (1982) formulation (Perry and Marusic eq.9) of the Coles wake function for
% managing mean boundary layer profiles.
%
% Syntax:  
%       [r13A, r13B] = getReynoldsStress13(eta, Pi, S, zeta, beta)
%       Determines reynolds stress 1 - 3: mean(u1'u3')/Utau^2 resulting from
%       type A and B eddies. Total Reynolds Stress is r13 = r13A + r13B.
%
% Inputs:
%
%       eta         [nZ x 1]    Nondimensional height z/deltac above the wall at
%                               which Reynolds Stresses are required. Eta values
%                               must ascend but not necessarily be monotonic.
%       
%   	Pi          [1 x 1]     Coles wake parameter Pi
%
%       S           [1 x 1]     Ratio between free stream and friction velocity
%                               S = U1/Utau
%
%       beta        [1 x 1]     The Clauser parameter, representing
%                               acceleration/decelaration of the boundary layer
%
%       zeta        [1 x 1]     Represents a scaled streamwise derivative of Pi
%
% Outputs:
%
%   	r13A        [nZ x 1]    Reynolds Stress \overline{u_1'u_3')/Utau^2
%                               resulting from Type A eddies.
%
%       r13B        [nZ x 1]    As r13A but representing the contribution from
%                               Type B eddies
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
%   [4] For the given wall formulation, can f be calculated more efficiently?
%       See schlichting and gersten p. 593
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
% Revisions:        30 July 2014        Created.
%                   18 April 2015       Properly documented with OAS header.
%                                       Changed PiParam variable name to Pi for
%                                       consistency with other functions in the
%                                       FCS.
%                   19 April 2015       Updated to use the modified Lewkowicz
%                                       formulation for Reynolds Stress, eq. 51
%                                       Perry and Marusic, rather than the
%                                       previously used logarithmic version (eq
%                                       53).
%
% Copyright (c) 2014-2015 Ocean Array Systems, All Rights Reserved.


% g1 = f2/S
% g2 = -f3/(C1*S)
% zeta = S*delta_c*(dPi/dx)
% C1 see eq 14

% define vk constant
kappa = 0.41;

% get f between eta = 0 and eta = 1 (bounds for C1 integration)
f      = getf(eta, Pi, kappa,S);
dPi = 0.01*Pi;
fplus  = getf(eta, Pi+dPi, kappa, S);
fminus = getf(eta, Pi-dPi, kappa, S);

% get C1
C1 = trapz(eta,f);



% Do the integrations for ei with central differencing for numerical
% differentiation of e coefficients 5-7
e1      = cumtrapz(eta,f);     
e1plus  = cumtrapz(eta,fplus);
e1minus = cumtrapz(eta,fminus);
e2      = cumtrapz(eta,f.^2);    
e2plus  = cumtrapz(eta,fplus.^2);
e2minus = cumtrapz(eta,fminus.^2);
e3 = f.*e1;                        
e4 = eta.*f;                       
e5 = (e1plus - e1minus)/(2*dPi);    
e6 = (e2plus - e2minus)/(2*dPi);
e7 = f.*e5;

% get A coefficients from equations A2a-d
A1 = e2 - e3 + S*(e4 - e1);
A2 = -2*e2 + e3 + S*e1;
A3 = e6 - e7 - S*e5;
A4 = 2*e2 - e3 + S*e4 - S*3*e1;

% B coefficients are simply evaluated at eta = 1 (eqn A6)
B1 = A1(end);
B2 = A2(end);
B3 = A3(end);
B4 = A4(end);

% E1 from (eqn A4)
E1 = 1/(kappa*S + 1);

% N from (eqn A5) using central differencing as before
wcminus = colesWake(1, Pi-dPi);
wcplus =  colesWake(1, Pi+dPi);
N = colesWake(1, Pi) + Pi*(wcplus-wcminus)/(2*dPi);

% Do the weird adding up thing for fi
f1 = 1 - A1./(B1 + E1.*B2)...
    - E1.*A2./(B1 + E1.*B2);

f2 = (E1.*N.*A2.*B1 ...
    + A3.*B1 - E1.*N.*A1.*B2 ...
    + E1.*A3.*B2 - A1.*B3 ...
    - E1.*A2.*B3) ./ (B1 + E1.*B2);

f3 = (E1.*A2.*B1...
    + A4.*B1 - E1.*A1.*B2...
    + E1.*A4.*B2 - A1.*B4...
    - E1.*A2.*B4) ./ (B1 + E1.*B2);

% convert f2 and f3 into g1 and g2 ready for application of eq. 15
g1 = f2/S;
g2 = -f3/(C1*S);

% Top 3 boxes of figure 18
minusReynStress = f1 + g1*zeta + g2*beta;

% Get the component due to equilibrium sink flow (OLD version - see P&M eqns
% 51,53)
% minusReynStressA = ones(size(eta)) - eta + eta.*log(eta);

% Lewkowicz 1982 shear stress for equilibrium sink flow, Perry and Marusic eqn. 51
minusReynStressA = 1 - (60/59)*eta - (20/59)*eta.^3 + (45/59)*eta.^4 - (24/59).*eta.^5 + (60/59)*eta.*log(eta);

minusReynStressA(isinf(minusReynStressA)) = 1; % Handle the log(0) singularity;

minusReynStressB = minusReynStress - minusReynStressA;

r13A = -1*minusReynStressA;

r13B = -1*minusReynStressB;

end


function [f] = getf(eta, Pi, kappa, S)

f = (-1/kappa)*log(eta) + (Pi/kappa)*colesWake(1,Pi)*ones(size(eta)) - (Pi/kappa)*colesWake(eta,Pi);
f(isinf(f)) = S; % from eqs 2 and 7 P&M part 1 1995

end 


function wc = colesWake(eta, Pi)

wc = 2*eta.^2.*(3-2*eta) - (1/Pi).*eta.^2.*(1-eta).*(1-2*eta);
end

