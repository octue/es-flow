function [] = basicGantt()
%BASICGANTT Summary of this function goes here
%   Detailed explanation goes here
 
% Number of elements in the discretised PDFS
nElems = 200;

% Dependencies
activity(1).predecessorsFinish = [];
activity(2).predecessorsFinish = [1];
activity(3).predecessorsFinish = [2];

% Properties (time in days)
activity(1).time_estimated = 48;
activity(2).time_estimated = 32;
activity(3).time_estimated = 18;

% Probability Density functions
activity(1).pdf = [10 0; 48 1; 62 0];
activity(2).pdf = [22 0; 32 1; 42 0];
activity(3).pdf = [10 0; 18 1; 50 0];
for actCtr = 1:numel(activity)
    pdf = activity(actCtr).pdf;
    pdf(:,2) = normalisePDF(pdf(:,1),pdf(:,2));
    pdf = rediscretisePDF(pdf,nElems);
end




end % End MAIN function



