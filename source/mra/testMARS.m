% Load a time series of ADCP data from about mid height in the water column
a = load('/Users/thc29/Documents/MATLAB/OAS/ADCP Data/EMEC/data.mat');
bin = 13;
ux = a.uxProf(:,bin);
uy = a.uyProf(:,bin);
uz = a.uzProf(:,bin);
disp(['Proportion of nan vectors in this bin: ' num2str(sum(isnan(ux))*100/numel(ux)) '%'])

% Clean up data to remove nans (linear interpolation where nans exist)
nanInds = isnan(ux) | isnan(uy) | isnan(uz);
indSet = 1:numel(ux);
truncIndSet = indSet(~nanInds);
ux = interp1(indSet(~nanInds), ux(~nanInds), indSet, 'linear', 0);
uy = interp1(indSet(~nanInds), uy(~nanInds), indSet, 'linear', 0);
uz = interp1(indSet(~nanInds), uz(~nanInds), indSet, 'linear', 0);
disp(['Proportion of nan vectors after cleaning: ' num2str(sum(isnan(ux))*100/numel(ux)) '%'])

% Clear up memory
clearvars a uxProf uyProf uzProf

% Trim number of samples
nSamp = 200000;
ux = ux(1:nSamp);
uy = uy(1:nSamp);
uz = uz(1:nSamp);