function [pdfi] = rediscretisePDF(x, prob, nElems)
%REDISCRETISEPDF Discretises a PDF with input number of samples using linear
%interpolation

xi = linspace(x(1),x(end),nElems)';
pdfi = interp1(x(:),prob(:),xi,'linear');
pdfi = [xi pdfi];

end

