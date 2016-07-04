

%% Load the test case and get the eddy spectra

load('eddyVelocities.mat')

i = 1;
[g, J, k1z, lambda] = getEddySpectra(types{i}, eddy(i).X,  eddy(i).Y,  eddy(i).Z,  eddy(i).U,  eddy(i).V,  eddy(i).W);







