function getEddySignatures(dir)
%GETEDDYSIGNATURE Computes all eddy signatures and saves them to signature_*.mat
% files in the input directory (default pwd)

%#ok<*NASGU>

if nargin == 0
    dir = pwd;
end

types = {'A', 'B1', 'B2', 'B3', 'B4'};

for i = 1:5
    
    % Get the type string
    type = types{i};
    
    % Get the eddy intensity functions
    [J, lambda, X, Y, Z, U, V, W, I] = getEddyIntensity(type); %#ok<ASGLU>

    % Get the spectra
    [g] = getEddySpectra(type, J, lambda, X, Y, Z, U, V, W);

    % Get the deficit
    [h] = getEddyDeficit(type, lambda);

    % Save spatial arrays to file in a more compact form to save memory:
    domain_extents = [X(1, 1, 1), X(1, end, 1)
                      Y(1, 1, 1), Y(end, 1, 1)
                      Z(1, 1, 1), Z(1, 1, end)];
    domain_spacing = [(X(1,2,1) - X(1,1,1))
                      (Y(2,1,1) - Y(1,1,1))
                      (Z(1,1,2) - Z(1,1,1))]';
    
    % QA check that no divide-by-0s slipped in
    if any(isnan(J(:))) || any(isinf(J(:)))
        throw('nan alert in J!!!')
    end
    if any(isnan(lambda(:))) || any(isinf(lambda(:)))
        throw('nan alert in lambda!!!')
    end
    if any(isnan(g(:))) || any(isinf(g(:)))
        throw('nan alert in g!!!')
    end
    % Save the file, overwriting any previous
    % TODO optionally save U,V,W,X,Y,Z
    file_name = fullfile(dir, ['signatures_' type '.mat']);
    save(file_name, 'J', 'type', 'lambda', 'domain_extents', 'domain_spacing', 'I', 'g', 'h', '-v7.3')
              
end

