function [ademResult] = adem(deltac, U1, Pi, S, zeta, beta)
%ADEM Runs the Attached-Detached Eddy Model for given boundary layer parameters
% Uses the using the Lewkowicz (1982) formulation (Perry and Marusic eq.9) of
% the Coles wake function to determine U(z), spectra and Reynolds Stresses from
% input parameters Pi, S, deltac, U1, zeta, beta.
%
% Syntax:  
%       [ademResult] = adem(z, deltac, U1, Pi, S, zeta, beta)
%       Determines all boundary layer parameters using the approach described in
%       Ref. 1.
%
% Inputs:
%       z           [nZ x 1]    Height in m above the wall at which output
%                               profiles are required
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
%       beta        [1 x 1]     The Clauser parameter, representing
%                               acceleration/decelaration of the boundary layer
%
%       zeta        [1 x 1]     Represents a scaled streamwise derivative of Pi
%
% Outputs:
%
%       adem        structure   Contains the following fields:
%
%           .Ux         [nZ x 1]    Streamwise velocity in m/s at points
%                                   corresponding to the heights in z
%
%           .z          [nZ x 1]    Heights above the wall
%
%           .Pi         [1 x 1]     As input
%           .S          [1 x 1]     As input
%           .deltac     [1 x 1]     As input
%           .U1         [1 x 1]     As input
%           .beta       [1 x 1]     As input
%           .zeta       [1 x 1]     As input
%
% Example:
%   
%   Replicate the Reynolds Stresses and Spectra from Skare and Krogstad 1994
%   [ademResult] = adem(1, 1, 6.85, 59.4, 0, 19)
%       
% References:
%
%   [1] Perry AE and Marusic I (1995) A wall-wake model for turbulent boundary
%       layers. Part 1. Extension of the attached eddy hypothesis J Fluid Mech
%       vol 298 pp 361-388
%
% Future Improvements:  
%
%   [1] See future improvements from all the other functions!
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
% Revisions:                19 April 2015       Created
%
% Copyright (c) 2014-2015 Ocean Array Systems, All Rights Reserved.


% Establish eddy intensity and spectral functions
eddyFile = 'eddySignatures.mat';
path = fileparts(which('adem'));
if exist(fullfile(path,eddyFile),'file')~=2
    [gA,  JA, k1z, lambda]  = getEddySpectra('A');
    [gB1, JB1] = getEddySpectra('B1');
    [gB2, JB2] = getEddySpectra('B2');
    [gB3, JB3] = getEddySpectra('B3');
    [gB4, JB4] = getEddySpectra('B4');
else
    s = load(eddyFile);
    k1z = s.eddy(1).k1z;
    lambda = s.eddy(1).lambda(:);
    gA = s.eddy(1).g;
    JA = s.eddy(1).J;
    gB1 = s.eddy(2).g;
    JB1 = s.eddy(2).J;
    gB2 = s.eddy(3).g;
    JB2 = s.eddy(3).J;
    gB3 = s.eddy(4).g;
    JB3 = s.eddy(4).J;
    gB4 = s.eddy(5).g;
    JB4 = s.eddy(5).J;
    clearvars s
end

% Ensemble average the Type B eddies
gB = (gB1+gB2+gB3+gB4)/4;
JB = (JB1+JB2+JB3+JB4)/4;

% Deconvolve for T^2w
[T2wA, T2wB, lambdaE, residualA, residualB] = getTw(Pi, S, zeta, beta, JA(:,3), JB(:,3), lambda);

% Get the mean profile at the same vertical points
eta = 1./exp(lambdaE);
z = eta*deltac;
[Ux] = getMeanProfile(Pi, S, deltac, U1, z);

% Determine Reynolds Stresses by convolution
[R, RA, RB] = getReynoldsStresses(T2wA, T2wB, JA, JB);

% Determine Spectra by convolution
[Psi, PsiA, PsiB] = getSpectra(T2wA, T2wB, gA, gB);

% Assemble results into output
ademResult.z            = z;
ademResult.eta          = eta;
ademResult.lambdaE      = lambdaE;
ademResult.deltac       = deltac;
ademResult.U1           = U1;
ademResult.Pi           = Pi;
ademResult.S            = S;
ademResult.zeta         = zeta;
ademResult.beta       	= beta;
ademResult.Ux           = Ux;
ademResult.eddyTypes    = {'A','B1','B2','B3','B4'};
ademResult.T2wA         = T2wA;
ademResult.T2wB         = T2wB;
ademResult.residualA    = residualA;
ademResult.residualB    = residualB;
ademResult.R            = R;
ademResult.RA           = RA;
ademResult.RB           = RB;
ademResult.Psi          = Psi;
ademResult.PsiA         = PsiA;
ademResult.PsiB         = PsiB;




end

