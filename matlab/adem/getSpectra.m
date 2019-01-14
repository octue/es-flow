function [Psi, PsiA, PsiB] = getSpectra(T2wA, T2wB, gA, gB, U1, S)
%GETSPECTRA Gets spectra 
%
% Syntax:  
%       [Psi, PsiA, PsiB] = getSpectra(T2wA, T2wB, gA, gB, U1, S)
%       Returns premultiplied boundary layer spectra according to Perry and
%       Marusic eq. 43. Divide these by k1z to get the power density spectrum 
%       Psi.
%
% Inputs:
%
%
% Outputs:       
%       
%   	Psi         [nZ x nk]       Total spectra as a function of height off
%                                   the wall and wavenumber k1z. 
%                                   Psi = PsiA + PsiB
%
%   	PsiA        [nZ x nk]       Contribution of Type A eddies to total
%                                   spectra as a function of height off 
%                                   the wall and wavenumber k1z.
%
%   	PsiA        [nZ x nk]       Contribution of Type B eddies to total
%                                   spectra as a function of height off 
%                                   the wall and wavenumber k1z.
%
% 
% References:
%
%   [1] Perry AE and Marusic I (1995) A wall-wake model for turbulent boundary
%       layers. Part 1. Extension of the attached eddy hypothesis J Fluid Mech
%       vol 298 pp 361-388
%
% Future Improvements:
%
% Author:                   T. Clark
% Work address:             Octue
%                           Hauser Forum
%                           3 Charles Babbage Road
%                           Cambridge
%                           CB3 0GT
% Email:                    tom@octue.com
% Website:                  www.octue.com
%
% Copyright (c) 2014-2018 Octue, All Rights Reserved.

% For each of the auto/cross spectra terms...
for j = 1:size(gA,3)
    % For each of the different heights...
    for i = 1:size(gA,2)
        PsiA(:,i,j) = conv(gA(:,i,j)', T2wA); %#ok<AGROW>
        PsiB(:,i,j) = conv(gB(:,i,j)', T2wB); %#ok<AGROW>
    end
end

% Trim the zero padded edges from the full convolution
PsiA = PsiA(2:end-1, :, :);
PsiB = PsiB(2:end-1, :, :);

% Remove the UTau^2 from eq. 43
Utau = U1/S;
PsiA = PsiA.*Utau.^2;
PsiB = PsiB.*Utau.^2;

% Summate from the different components
Psi = PsiA + PsiB;
