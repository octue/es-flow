function [filtered, residual] = extractCoherent(signal, varargin)
%EXTRACTCOHERENT Remove the first dyadic level (decorrelated noise)
%   Uses biorthogonal (2:2) wavelet decomposition to limit scales to those which
%   can be expressed to second order smoothness with the sampling rate
%
% Syntax:
%       [filtered] = extractCoherent(signal) Returns approximation to the signal
%
%       [filtered, residual] = extractCoherent(signal) Returns residual
%       components too
%
%       [filtered, residual] = extractCoherent(signal, N) Removes up the the Nth
%       dyadic level
%
%
% Inputs:
%
%       signal          [n x 3]     n point velocity vectors [u v w]
%
%       N               [1 x 1]     integer, default 1 (recommended)
%
% Outputs:
%
%       filtered        [n x 3]     n point filtered approximation to the input
%
%       residual        [n x 3]     n point residual
%
% See Also: WAVEMENU
%
% Future Improvements:
%
%   [1] Addition of lengthscales based on signal content or dyadic level to
%       remove
%
% References:
%
% Author:                   T. H. Clark
% Work address:             Ocean Array Systems Ltd
%                           Hauser Forum
%                           3 Charles Babbage Road
%                           Cambridge
%                           CB3 0GT
% Email:                    tom.clark@oceanarraysystems.com
% Website:                  www.oceanarraysystems.com
%
% Copyright (c) 2015-6 Ocean Array Systems, All Rights Reserved.

% Allow for multiple levels to be taken
if nargin == 2
    N = varargin{1};
else
    % Use only one level of decomposition
    N = 1;
end

% Decompose with a Biorthogonal wavelet, second order
[Cu, Lu] = wavedec(signal(:,1), N, 'bior2.2');
[Cv, Lv] = wavedec(signal(:,2), N, 'bior2.2');
[Cw, Lw] = wavedec(signal(:,3), N, 'bior2.2');

% Split into approximation and detail coefficients
Au = appcoef(Cu, Lu, 'bior2.2', N);
Av = appcoef(Cv, Lv, 'bior2.2', N);
Aw = appcoef(Cw, Lw, 'bior2.2', N);
Bu = detcoef(Cu, Lu, 'bior2.2', N);
Bv = detcoef(Cv, Lv, 'bior2.2', N);
Bw = detcoef(Cw, Lw, 'bior2.2', N);

uFilt = waverec([Au; zeros(size(Bu{1}))], Lu, 'bior2.2');
vFilt = waverec([Av; zeros(size(Bv{1}))], Lv, 'bior2.2');
wFilt = waverec([Aw; zeros(size(Bw{1}))], Lw, 'bior2.2');

uRes = waverec([zeros(size(Au)); Bu{1}], Lu, 'bior2.2');
vRes = waverec([zeros(size(Av)); Bv{1}], Lv, 'bior2.2');
wRes = waverec([zeros(size(Aw)); Bw{1}], Lw, 'bior2.2');

% Return filtered and residual components - we just get rid of the whole part
filtered = [uFilt   vFilt   wFilt];
residual = [uRes    vRes    wRes];

end

