function [Jz, J, lambda, X, Y, Z, U, V, W] = getEddySpectra_v(varargin)
%GETEDDYSPECTRA Returns eddy spectral functions gij, with i=1:3, j=1:3
% These are the signatures (in terms of turbulent fluctuations) that an
% individual structure will contribute to a boundary layer. Also calculates
% intensity functions J.
%
% Syntax:
%       [g, J, k1z, lambda] = getEddySpectra('type')
%       Returns eddy spectral and intensity functions gij for the eddy type
%       described in the input string, presently 'A', 'B1', 'B2', 'B3', 'B4'
%       which are the eddies described in Ref 1.
%
%       [..., X, Y, Z, U, V, W] = getEddySpectra('type')
%       Also returns the velocity signature of the eddy on a regularised grids
%       (as produced by meshgrid).
%
%       [...] = getEddySpectra([J, lambda, X, Y, Z, U, V, W])
%       Performs just the eddy spectral calculations without recomputing
%       intensity functions.
%
% Inputs:
%
%       type            string      The eddy type for which to calculate J.
%                                   Current accepted values are  'A', 'B1',
%                                   'B2', 'B3', 'B4', which correspond with the
%                                   eddy types in Ref [1].
%
% Outputs:
%       g               [nZ x nk x 6]
%                                   g(:,:,page) contains the eddy spectral
%                                   function gij which is [nZ x nk] in size.
%                                   Array pages contain gij ordered as:
%                                       [g11 g12 g13 g22 g23 g33];
%
%       J               [nZ x 6]    Jij contains the eddy intensity
%                                   function Jij(lambda), which is [nZ x 1] in
%                                   size. Note that the lower diagonal is not
%                                   included due to symmetry of the Reynolds
%                                   Stress Tensor. Array columns contain Jij,
%                                   ordered as:
%                                       [J11 J12 J13 J22 J23 J33];
%
%       k1z             [nZ x nK]   Wavenumber used to premultiply the power
%                                   spectra
%
%       lambda          [nZ x 1]    The remapped wall coordinate at which J is
%                                   given. Varies from 0 (at the edge of the BL)
%                                   to a point close to the wall.
%
%       X,Y,Z           [nY, nX, nZ]
%                                   Regular arrays as produced by meshgrid
%                                   containing the volume coordinate points at
%                                   which the velocity signaturee was computed.
%
%       U,V,W           [nY, nX, nZ]
%                                   Velocities in m/s induced by a single eddy
%                                   of scale z/delta = 1 and of strength gamma =
%                                   1, together with its image in the wall (z =
%                                   0), computed on the regular grid in X,Y,Z.
%
%
% See Also: GETEDDYINTENSITY.M
%
% Future Improvements:
%
%   [1] Computation could be done using the eddy spectral function gij (eq.44
%       Perry and Marusic 1995) which is quicker but more complicated than the
%       direct integration of eqn 35 that is presently implemented.
%
%   [2] Presently the intensity functions are computed each time this is called.
%       There should be a cached lookup for particular eddy types (i.e. a .mat
%       file storing all Jij).
%
%   [3] Check why we're discretising the line filaments into small pieces -
%       probably a hangover from the old biot savart code. Eliminating (setting
%       nEl = 1) could save a lot of compute time right there (factor 51
%       improvement).
%
%   [4] Possible additional modification to take into account a free surface
%       image.
%
% References:
%
%   [1] Perry AE and Marusic I (1995) A wall-wake model for turbulent boundary
%       layers. Part 1. Extension of the attached eddy hypothesis J Fluid Mech
%       vol 298 pp 361-388
%
% Author:                   T. H. Clark (Vik - modified on 18th May)
% Work address:             Ocean Array Systems Ltd
%                           Hauser Forum
%                           3 Charles Babbage Road
%                           Cambridge
%                           CB3 0GT
% Email:                    tom.clark@oceanarraysystems.com
% Website:                  www.oceanarraysystems.com
%
% Revision History:        	19 April 2015       Created from the deprecated getJ
%                                               Added the string type specifier
%                                               and reworked for new biot savart
%                                               interface.
%                           23 April 2014       Added volumetric grid and
%                                               velocity signature to the output
%                           24 April 2014       Added documentation for the full
%                                               input of intensity functions
%
% Copyright (c) 2014-2015 Ocean Array Systems, All Rights Reserved.

if numel(varargin) == 1
    % Use the getEddyIntensity function to get velocity arrays
    [J, lambda, X, Y, Z, U, V, W] = getEddyIntensity(varargin{1});
else
    % Load already computed velocity arrays from the input (allows for quicker
    % modification and debug of this function without recomputation of all the
    % eddy signatures)
    [J, lambda, X, Y, Z, U, V, W] = deal(varargin{:}); clearvars varargin
end

% Interpolate UVW to the same locations in lambda that we have the J functions
% oldZVec = Z(1,1,:);
% newZVec = (1./exp(lambda));
% newU = zeros(size(U,1), size(U,2), numel(newZVec));
% newV = zeros(size(U,1), size(U,2), numel(newZVec));
% newW = zeros(size(U,1), size(U,2), numel(newZVec));
% dispstat('Determining eddy spectral function...','init','timestamp')
% for i = 1:size(U,1);
%     dispstat(['Determining eddy spectral function, profile ' num2str(i) ' of ' num2str(size(U,1)) '.'],'timestamp')
%     for j = 1:size(U,2)        
%         uvec = U(i,j,:);
%         vvec = V(i,j,:);
%         wvec = W(i,j,:);
%         newU(i,j,:) = interp1(oldZVec(:), uvec(:), newZVec(:), 'pchip');
%         newV(i,j,:) = interp1(oldZVec(:), vvec(:), newZVec(:), 'pchip');
%         newW(i,j,:) = interp1(oldZVec(:), wvec(:), newZVec(:), 'pchip');
%     end
% end

% save
% load

%% Take the fourier transform of each component along the streamwise direction
% (dimension 2 in this case)
% Fi = fft(newU, [], 2); %clearvars newU
% Fj = fft(newV, [], 2); %clearvars newV
% Fk = fft(newW, [], 2); %clearvars newW
% 
% % Kill the redundant half of the FFT?
% Fi(:,end/2:end,:) = 0;
% Fj(:,end/2:end,:) = 0;
% Fk(:,end/2:end,:) = 0;
% 
% % Scale by N
% Fi = Fi./size(Fi,2);
% Fj = Fj./size(Fi,2);
% Fk = Fk./size(Fi,2);

k1 = 0:0.1:10;
Nk = length(k1);
Fio = zeros(size(X,1),Nk,size(X,3));
Fjo = zeros(size(X,1),Nk,size(X,3));
Fko = zeros(size(X,1),Nk,size(X,3));

for ii = 1:Nk
    
    Fio(:,ii,:) = trapz(X(1,:,1),U.*exp(-2i*pi*k1(ii)*X),2);
    Fjo(:,ii,:) = trapz(X(1,:,1),V.*exp(-2i*pi*k1(ii)*X),2);
    Fko(:,ii,:) = trapz(X(1,:,1),W.*exp(-2i*pi*k1(ii)*X),2);
    
end

J11z = trapz(k1,trapz(Y(:,1,1),real((conj(Fio).*Fio)),1),2);
J12z = trapz(k1,trapz(Y(:,1,1),real((conj(Fio).*Fjo)),1),2);
J13z = trapz(k1,trapz(Y(:,1,1),real((conj(Fio).*Fko)),1),2);
J22z = trapz(k1,trapz(Y(:,1,1),real((conj(Fjo).*Fjo)),1),2);
J23z = trapz(k1,trapz(Y(:,1,1),real((conj(Fjo).*Fko)),1),2);
J33z = trapz(k1,trapz(Y(:,1,1),real((conj(Fko).*Fko)),1),2);
%%

Jz   = zeros(6,size(X,3));

for ii = 1:size(X,3)

Jz(:,ii) = [J11z(1,1,ii); J12z(1,1,ii); J13z(1,1,ii); J22z(1,1,ii); J23z(1,1,ii); J33z(1,1,ii)];

end
%%
% Get the k1delta
% Dx = X(1,2,1) - X(1,1,1);
% N = size(X,2);
% k1delta = (0:N-1)*2*pi/Dx;
% 
% % Integrate over the y direction for Gij (eqn 40)
% zVec = reshape(newZVec,[1 1 numel(newZVec)]);
% deltaOnz = repmat(1./zVec, [1 size(Z,2) 1]);
% G11 = deltaOnz .* trapz(Y(:,1,1), real(conj(Fi).*Fi), 1);
% G12 = deltaOnz .* trapz(Y(:,1,1), real(conj(Fi).*Fj), 1);
% G13 = deltaOnz .* trapz(Y(:,1,1), real(conj(Fi).*Fk), 1); clearvars Fi
% G22 = deltaOnz .* trapz(Y(:,1,1), real(conj(Fj).*Fj), 1);
% G23 = deltaOnz .* trapz(Y(:,1,1), real(conj(Fj).*Fk), 1); clearvars Fj
% G33 = deltaOnz .* trapz(Y(:,1,1), real(conj(Fk).*Fk), 1); clearvars Fk
% 
% 
% if true
%     raiseFigure('G11'); contourf(permute(G11,[2 3 1]),30); axis equal; colorbar; title('G11(1,i,j)'); xlabel('i'); ylabel('j');
%     raiseFigure('G12'); contourf(permute(G12,[2 3 1]),30); axis equal; colorbar; title('G12(1,i,j)'); xlabel('i'); ylabel('j');
%     raiseFigure('G13'); contourf(permute(G13,[2 3 1]),30); axis equal; colorbar; title('G13(1,i,j)'); xlabel('i'); ylabel('j');
%     raiseFigure('G22'); contourf(permute(G22,[2 3 1]),30); axis equal; colorbar; title('G22(1,i,j)'); xlabel('i'); ylabel('j');
%     raiseFigure('G23'); contourf(permute(G23,[2 3 1]),30); axis equal; colorbar; title('G23(1,i,j)'); xlabel('i'); ylabel('j');
%     raiseFigure('G33'); contourf(permute(G33,[2 3 1]),30); axis equal; colorbar; title('G33(1,i,j)'); xlabel('i'); ylabel('j');
% end
%     
% % Determine k1z (k1z = k1delta * z/delta)
% k1z = repmat(k1delta,[1 1 size(deltaOnz,3)])./deltaOnz;
% 
% % Multiply for gij
% g11 = k1z.*G11;
% g12 = k1z.*G12;
% g13 = k1z.*G13;
% g22 = k1z.*G22;
% g23 = k1z.*G23;
% g33 = k1z.*G33;
% 
% if true
%     raiseFigure('g11'); contourf(permute(g11,[2 3 1]),30); axis equal; colorbar; title('g11(1,i,j)'); xlabel('i'); ylabel('j');
%     raiseFigure('g12'); contourf(permute(g12,[2 3 1]),30); axis equal; colorbar; title('g12(1,i,j)'); xlabel('i'); ylabel('j');
%     raiseFigure('g13'); contourf(permute(g13,[2 3 1]),30); axis equal; colorbar; title('g13(1,i,j)'); xlabel('i'); ylabel('j');
%     raiseFigure('g22'); contourf(permute(g22,[2 3 1]),30); axis equal; colorbar; title('g22(1,i,j)'); xlabel('i'); ylabel('j');
%     raiseFigure('g23'); contourf(permute(g23,[2 3 1]),30); axis equal; colorbar; title('g23(1,i,j)'); xlabel('i'); ylabel('j');
%     raiseFigure('g33'); contourf(permute(g33,[2 3 1]),30); axis equal; colorbar; title('g33(1,i,j)'); xlabel('i'); ylabel('j');
% 
%     %raiseFigure('eddy ribbons'); clf; streamribbon(X,Y,Z,U,V,W,rand(10,1),(rand(10,1)-0.5)*2,rand(10,1)*1.5); shading interp; view(3); camlight; lighting gouraud;
% end
% 
% % These are [1 x nk x nZ] arrays and we should probably put them in a reusable
% % form. Swap dimensions to [nZ x nk]:
% order = [3 2 1];
% k1z = permute(k1z,order);
% g11 = permute(g11,order);
% g12 = permute(g12,order);
% g13 = permute(g13,order);
% g22 = permute(g22,order);
% g23 = permute(g23,order);
% g33 = permute(g33,order);
% 
% % Output to an array concatenated in the third dimension
% g = cat(3, g11, g12, g13, g22, g23, g33);
% 
% % We should now be able to retrieve the Jij functions by integrating g wrt
% % alphaz as per eqn. 44. Annoyingly there is a -Inf value at the zero
% % wavenumbers in alphaz.
% alphaz = log(k1z);
% alphaz(:,1) = alphaz(:,2);
% da = diff(alphaz,1,2);
% checkJ11 = trapz(da.*g11(:,2:end),2);
% checkJ12 = trapz(da.*g12(:,2:end),2);
% checkJ13 = trapz(da.*g13(:,2:end),2);
% checkJ22 = trapz(da.*g22(:,2:end),2);
% checkJ23 = trapz(da.*g23(:,2:end),2);
% checkJ33 = trapz(da.*g33(:,2:end),2);
% checkJ = [checkJ11 checkJ12 checkJ13 checkJ22 checkJ23 checkJ33];
% 
% raiseFigure('Check J between direct calculation and spectral'); 
% clf; 
% lh = plot(lambda, checkJ,'--'); 
% hold on; 
% lh2 = plot(lambda,J,'-');
% for i = 1:numel(lh2)
%     set(lh2(i),'Color',get(lh(i),'Color'))
% end
% hold on;
% 
% legend({'J11 direct';
%     'J12 direct';
%     'J13 direct';
%     'J22 direct';
%     'J23 direct';
%     'J33 direct';
%     'J11 spectral';
%     'J12 spectral';
%     'J13 spectral';
%     'J22 spectral';
%     'J23 spectral';
%     'J33 spectral'});