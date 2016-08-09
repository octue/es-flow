function [pdf] = linearPDF(x,a,b,c,d)
%LINEARPDF linearly varying probability density function
%  PDF = linearPDF(X, A, B, C, D) computes the probability density
%  function (PDF) at X of the of a linear probability function distributed
%  between A and B with end probabilities C and D. Note that the function is
%  normalised to give an cumulative probability of 1 over the range, 0 outside.



%% ERROR CHECKING
if nargin ~= 5
    error('MATLAB:monteCarlo:linearPDF','linearPDF: Incorrect number of arguments');
end

if ~isscalar(a) || ~isscalar(b) || ~isscalar(c) || ~isscalar(d)
    error('MATLAB:monteCarlo:linearPDF','linearPDF: Arguments A and B and C and D must be scalar parameters (location1, location2, probability1, probability2)');
end


%% CREATE DISTRIBUTION
lowEnd = min([a b]);
if lowEnd == a
    pdf = 



