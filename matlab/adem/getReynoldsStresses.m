function [R, RA, RB] = getReynoldsStresses(T2wA, T2wB, JA, JB)
%GETREYNSTRESSES Outputs all components of Reynolds Stresses


% For each of the different Reynolds Stress terms... 
for i = 1:size(JA,2)
    
    RA(:,i) = conv(JA(3:end,i), T2wA); %#ok<AGROW>
    RB(:,i) = conv(JB(3:end,i), T2wB); %#ok<AGROW>
    
end

% Add to get full Reynolds Stresses
R = RA + RB;

% Plot Reynolds Stress distributions for checking purposes
raiseFigure('Reynolds Stress (A)')
clf
plot(RA)

raiseFigure('Reynolds Stress B')
clf
plot(RB)
