

%% Load the test case and get the eddy spectra

load('eddySignatures.mat')

i = 1;
[Jz, J, lambda, X, Y, Z, U, V, W] = getEddySpectra_vik(eddy(i).J, eddy(i).lambda, eddy(i).X,  eddy(i).Y,  eddy(i).Z,  eddy(i).U,  eddy(i).V,  eddy(i).W);

z = zeros(1,size(Z,3));
for jj = 1:size(Z,3)
    z(jj) = Z(1,1,jj);
end

lambda1 = log(1./z);

figure
plot(lambda,J(:,1),'k')
hold on
plot(lambda,J(:,2),'r')
plot(lambda,J(:,3),'b')
plot(lambda,J(:,4),'c')
plot(lambda,J(:,5),'g')
plot(lambda,J(:,6),'m')

figure
plot(lambda1,Jz(1,:),'o-k')
hold on
plot(lambda1,Jz(2,:),'o-r')
plot(lambda1,Jz(3,:),'o-b')
plot(lambda1,Jz(4,:),'o-c')
plot(lambda1,Jz(5,:),'o-g')
plot(lambda1,Jz(6,:),'o-m')














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





