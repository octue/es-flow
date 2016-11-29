

%% Load the test case and get the eddy spectra

load('eddySignatures.mat')

i = 1;
[g, ~, k1z] = getEddySpectra(eddy(i).J, eddy(i).lambda, eddy(i).X,  eddy(i).Y,  eddy(i).Z,  eddy(i).U,  eddy(i).V,  eddy(i).W);





















% %% Create reduced test case
% load('eddySignatures.mat')
% for i = 1:5
%     
% 	test(i).lambda = eddy(i).lambda(1:4:end);
%     test(i).J = eddy(i).J(1:4:end,:);
% 	test(i).h = eddy(i).h(1:4:end,:);
% 
% 	test(i).X = eddy(i).X;
% 	test(i).Y = eddy(i).Y;
% 	test(i).Z = eddy(i).Z;
% 	test(i).U = eddy(i).U;
% 	test(i).V = eddy(i).V;
% 	test(i).W = eddy(i).W;
% 	test(i).g = eddy(i).g;
% 	test(i).k1z = eddy(i).k1z;
% end
% 
% eddy = test;
% save('testSignatures.mat','eddy','types','-v7.3')
% 
% clear all





