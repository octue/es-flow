function [Psi, PsiA, PsiB] = getSpectra(T2wA, T2wB, gA, gB)
%GETSPECTRA Gets spectra 
%
% Syntax:  
%       [Psi, PsiA, PsiB] = getSpectra(T2wA, T2wB, gA, gB)
%       Returns premultiplied boundary layer spectra according to Perry and
%       Marusic eq. 43.
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
% References:
%
%   [1] Perry AE and Marusic I (1995) A wall-wake model for turbulent boundary
%       layers. Part 1. Extension of the attached eddy hypothesis J Fluid Mech
%       vol 298 pp 361-388
%
% Future Improvements:
%
% Author:                   T. Clark
% Work address:             Ocean Array Systems Ltd
%                           Hauser Forum
%                           3 Charles Babbage Road
%                           Cambridge
%                           CB3 0GT
% Email:                    tom.clark@oceanarraysystems.com
% Website:                  www.oceanarraysystems.com
%
% Created:                  19 April 2015       Created
%
% Copyright (c) 2014-2015 Ocean Array Systems, All Rights Reserved.


% Preallocate for speed
% PsiA = zeros(size(gA,2), size(gA,1)+numel(T2wA), size(gA,3));
% PsiA = zeros(size(gA,2), size(gA,1)+numel(T2wA), size(gA,3));

% For each of the auto/cross spectra terms...
for j = 1:size(gA,3)
    % For each of the different heights...
    for i = 1:size(gA,2)
        PsiA(:,i,j) = conv(gA(:,i,j)', T2wA);
        PsiB(:,i,j) = conv(gB(:,i,j)', T2wB);
    end
end

Psi = PsiA + PsiB;