
% Take the fourier transform of each component along the streamwise direction
% (dimension 2 in this case)
Fi = fft(newU, [], 2); %clearvars newU
Fj = fft(newV, [], 2); %clearvars newV
Fk = fft(newW, [], 2); %clearvars newW

% Kill the redundant half of the FFT?
Fi(:,end/2:end,:) = 0;
Fj(:,end/2:end,:) = 0;
Fk(:,end/2:end,:) = 0;

% Scale by N
Fi = Fi./size(Fi,2);
Fj = Fj./size(Fi,2);
Fk = Fk./size(Fi,2);

% Get the k1delta
Dx = X(1,2,1) - X(1,1,1);
N = size(X,2);
k1delta = (0:N-1)*2*pi/Dx;

% Integrate over the y direction for Gij (eqn 40)
zVec = reshape(newZVec,[1 1 numel(newZVec)]);
deltaOnz = repmat(1./zVec, [1 size(Z,2) 1]);
G11 = deltaOnz .* trapz(Y(:,1,1), real(conj(Fi).*Fi), 1);
G12 = deltaOnz .* trapz(Y(:,1,1), real(conj(Fi).*Fj), 1);
G13 = deltaOnz .* trapz(Y(:,1,1), real(conj(Fi).*Fk), 1); clearvars Fi
G22 = deltaOnz .* trapz(Y(:,1,1), real(conj(Fj).*Fj), 1);
G23 = deltaOnz .* trapz(Y(:,1,1), real(conj(Fj).*Fk), 1); clearvars Fj
G33 = deltaOnz .* trapz(Y(:,1,1), real(conj(Fk).*Fk), 1); clearvars Fk


if true
    raiseFigure('G11'); contourf(permute(G11,[2 3 1]),30); axis equal; colorbar; title('G11(1,i,j)'); xlabel('i'); ylabel('j');
    raiseFigure('G12'); contourf(permute(G12,[2 3 1]),30); axis equal; colorbar; title('G12(1,i,j)'); xlabel('i'); ylabel('j');
    raiseFigure('G13'); contourf(permute(G13,[2 3 1]),30); axis equal; colorbar; title('G13(1,i,j)'); xlabel('i'); ylabel('j');
    raiseFigure('G22'); contourf(permute(G22,[2 3 1]),30); axis equal; colorbar; title('G22(1,i,j)'); xlabel('i'); ylabel('j');
    raiseFigure('G23'); contourf(permute(G23,[2 3 1]),30); axis equal; colorbar; title('G23(1,i,j)'); xlabel('i'); ylabel('j');
    raiseFigure('G33'); contourf(permute(G33,[2 3 1]),30); axis equal; colorbar; title('G33(1,i,j)'); xlabel('i'); ylabel('j');
end
    
% Determine k1z (k1z = k1delta * z/delta)
k1z = repmat(k1delta,[1 1 size(deltaOnz,3)])./deltaOnz;

% Multiply for gij
g11 = k1z.*G11;
g12 = k1z.*G12;
g13 = k1z.*G13;
g22 = k1z.*G22;
g23 = k1z.*G23;
g33 = k1z.*G33;

if true
    raiseFigure('g11'); contourf(permute(g11,[2 3 1]),30); axis equal; colorbar; title('g11(1,i,j)'); xlabel('i'); ylabel('j');
    raiseFigure('g12'); contourf(permute(g12,[2 3 1]),30); axis equal; colorbar; title('g12(1,i,j)'); xlabel('i'); ylabel('j');
    raiseFigure('g13'); contourf(permute(g13,[2 3 1]),30); axis equal; colorbar; title('g13(1,i,j)'); xlabel('i'); ylabel('j');
    raiseFigure('g22'); contourf(permute(g22,[2 3 1]),30); axis equal; colorbar; title('g22(1,i,j)'); xlabel('i'); ylabel('j');
    raiseFigure('g23'); contourf(permute(g23,[2 3 1]),30); axis equal; colorbar; title('g23(1,i,j)'); xlabel('i'); ylabel('j');
    raiseFigure('g33'); contourf(permute(g33,[2 3 1]),30); axis equal; colorbar; title('g33(1,i,j)'); xlabel('i'); ylabel('j');

    %raiseFigure('eddy ribbons'); clf; streamribbon(X,Y,Z,U,V,W,rand(10,1),(rand(10,1)-0.5)*2,rand(10,1)*1.5); shading interp; view(3); camlight; lighting gouraud;
end

% These are [1 x nk x nZ] arrays and we should probably put them in a reusable
% form. Swap dimensions to [nZ x nk]:
order = [3 2 1];
k1z = permute(k1z,order);
g11 = permute(g11,order);
g12 = permute(g12,order);
g13 = permute(g13,order);
g22 = permute(g22,order);
g23 = permute(g23,order);
g33 = permute(g33,order);

% Output to an array concatenated in the third dimension
g = cat(3, g11, g12, g13, g22, g23, g33);

% We should now be able to retrieve the Jij functions by integrating g wrt
% alphaz as per eqn. 44. Annoyingly there is a -Inf value at the zero
% wavenumbers in alphaz.
alphaz = log(k1z);
alphaz(:,1) = alphaz(:,2);
da = diff(alphaz,1,2);
checkJ11 = trapz(da.*g11(:,2:end),2);
checkJ12 = trapz(da.*g12(:,2:end),2);
checkJ13 = trapz(da.*g13(:,2:end),2);
checkJ22 = trapz(da.*g22(:,2:end),2);
checkJ23 = trapz(da.*g23(:,2:end),2);
checkJ33 = trapz(da.*g33(:,2:end),2);
checkJ = [checkJ11 checkJ12 checkJ13 checkJ22 checkJ23 checkJ33];

raiseFigure('Check J between direct calculation and spectral'); 
clf; 
lh = plot(lambda, checkJ,'--'); 
hold on; 
lh2 = plot(lambda,J,'-');
for i = 1:numel(lh2)
    set(lh2(i),'Color',get(lh(i),'Color'))
end
hold on;

legend({'J11 direct';
    'J12 direct';
    'J13 direct';
    'J22 direct';
    'J23 direct';
    'J33 direct';
    'J11 spectral';
    'J12 spectral';
    'J13 spectral';
    'J22 spectral';
    'J23 spectral';
    'J33 spectral'});