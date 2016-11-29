function setupEddyDatabase()
%SETUPEDDYDATABASE Creates a .mat file containing eddy intensity functions and
%signatures, avoiding the need for repeated calculation.

types = {'A','B1','B2','B3','B4'};

if exist('eddyVelocities.mat','file') == 2
    s = load('eddyVelocities.mat');
    for iType = 1:numel(types)
        % Build the same structure containing eddy deficits, spectral functions
        % and intensity functions, but use the existing velocity arrays
        [eddy(iType).g, eddy(iType).J, eddy(iType).k1z, eddy(iType).lambda] = getEddySpectra(types{iType},...
            s.eddy(iType).X, s.eddy(iType).Y, s.eddy(iType).Z, ...
            s.eddy(iType).U, s.eddy(iType).V, s.eddy(iType).W);
        [eddy(iType).h] = getEddyDeficit(types{iType}, eddy(iType).lambda);
    
    end
else
    for iType = 1:numel(types)
        % Build a structure containing eddy deficits, spectral functions and
        % intensity functions along with the velocity signature arrays used to
        % calculate all of this stuff.
        [eddy(iType).g, eddy(iType).J, eddy(iType).k1z, eddy(iType).lambda,...
            eddy(iType).X, eddy(iType).Y, eddy(iType).Z, ...
            eddy(iType).U, eddy(iType).V, eddy(iType).W] = getEddySpectra(types{iType});
        [eddy(iType).h] = getEddyDeficit(types{iType}, eddy(iType).lambda);
    end
end

% Get destination folder
[path] = fileparts(which('adem'));

% Save the eddies and their types to a mat file in the toolbox directory
save(fullfile(path,'eddySignatures.mat'),'eddy','types','-v7.3');





