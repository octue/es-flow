function [T2wA, T2wB, lambdaE, residualJA, residualJB] = getTw(Pi, S, zeta, beta, J13A, J13B, lambda)
%GETTW Gets the T^2w distribution for calculation of Spectra and Stress terms
%
% Syntax:  
%       [T2wA, T2wB, lambdaE] = getT2w(Pi, S, zeta, beta, JA13, JB13)
%       Determines T^2w functions for Type aA and B eddies
%
%       [T2wA, T2wB, lambdaE, residualJA, residualJB] = getT2w(Pi, S, zeta, beta, JA13, JB13)
%       Also outputs residuals from the convolution to check
%
% Inputs:
%       
%   	Pi          [1 x 1]     Coles wake parameter Pi
%
%       S           [1 x 1]     Ratio between free stream and friction velocity
%                               S = U1/Utau
%
%       zeta        [1 x 1]     Represents a scaled streamwise derivative of Pi
%
%       beta        [1 x 1]     The Clauser parameter, representing
%                               acceleration/decelaration of the boundary layer
%
%       JA13        [nL x 1]    Contribution of a Type A eddy to the Reynolds
%                               Stress at points in lambda
%
%       JA13        [nL x 1]    Contribution of a Type A eddy to the Reynolds
%                               Stress at points in lambda
%
%       lambda      [nL x 1]    Logarithmic normalised wall coordinates Lambda
%                               at which J13A,B are provided.
%
% Outputs:
%
%       T           [ x ]       A measure of how U1/Utau varies with eddy scale;
%                               this is the strength PDF
%
%       w           [ x ]       A measure of how the PDF of eddy scales departs
%                               from a -1 power law; this is the scale PDF.
%
%       T2w         [ x ]       -1*(T.^2).*w
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
%   [5] Numerical quantity check on residuals of the deconvolution
%
%   [6] The shear distribution dUds/dEta is based on the lewkowicz 1982
%       formulation of the wake deficit function and needs to be changed along
%       with any modifications to the shear profile.
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
% Revisions:        unknown             Created.
%                   18 April 2015       Properly documented with OAS header.
%                                       Changed PiParam variable name to Pi for
%                                       consistency with other functions in the
%                                       FCS.
%                   19 April 2015       Altered to provide only T^2w fo rthe
%                                       time being
%
% Copyright (c) 2014-2015 Ocean Array Systems, All Rights Reserved.

% Method:
%   Get Shear and Reynolds Stress profiles for the present flow
%   Get Tw for A and B types by deconvolution of h from dUd*/dlambda
%   Get T^2w for A and B types by deconvolution of J13 from R13
%   Get T and w by element wise division
%   Perform the log conversion to get Q and D (the functions of strength and PDF
%   of scales)


%% GET SHEAR AND REYNOLDS STRESS PROFILES

% Define a range for lamdaE (lamdaE = ln(delta_c/z) = ln(1/eta))
lambdaE = 0:0.01:100;

% Re-express as eta and flip so that eta ascends (required for
% getReynoldsStress13)
eta = fliplr(exp(-1*lambdaE));
disp('eta 1')
disp(eta(1))
disp('eta end')
disp(eta(end))
% Get the Reynolds Stresses and flip back
[r13A, r13B] = getReynoldsStress13(eta, Pi, S, zeta, beta);
r13A = fliplr(r13A);
r13B = fliplr(r13B);


% %% GET Tw FOR A AND B EDDIES BY DECONVOLUTION
% 
% % Deconvolve out the A and B structure contributions to the mean shear
% [minusTwA, residualhA] = deconv(r13A, hA);
% [minusTwB, residualhB] = deconv(r13B, hB);



%% GET T^2w FOR A AND B EDDIES BY DECONVOLUTION

J13A = J13A(3:end);
J13B = J13B(3:end);


% Deconvolve out the A and B structure contributions to the Reynolds Stresses.
% NOTE: it's actually -1*T^2w that comes out.
disp('SIZES')
disp(size(r13A))
disp(size(J13A'))
[T2wA, residualJA] = deconv(r13A, J13A');
[T2wB, residualJB] = deconv(r13B, J13B');

raiseFigure('t2wA')
clf
subplot(1,3,1)
plot(lambdaE,r13A); hold on; plot(lambdaE, r13B)
legend({'R13A'; 'R13B'})
xlabel('\lambda_E')

subplot(1,3,2)
plot(J13A)
% plot(lambda(3:end),J13A)
hold on 
plot(J13B)
% plot(lambda(3:end),J13B)
legend({'J13A'; 'J13B'})
xlabel('\lambda')

subplot(1,3,3)
plot(T2wA)
hold on
plot(T2wB)
legend({'T^2\omegaA'; 'T^2\omegaB'})


%% CHECK CONVOLUTIONS
% Reconvolve to check we can go round without losing accuracy. TODO: Numerical
% quantity check. 
r13A_recon = conv(J13A',T2wA);
r13B_recon = conv(J13B',T2wB);

raiseFigure('getTw() Convolution Check')
clf
plot(lambdaE(:), [r13A(:), r13A_recon(:), r13B(:), r13B_recon(:)]);

legend({'R13A'; 'R13A reconstructed'; 'R13B'; 'R13B reconstructed'})




end

