function [psd, f] = spectrumData(fs, data)
%SPECTRUMDATA Returns One Sided PSD computed using direct FFT method for
%uniformly sampled data (i.e. NOT estimated as a periodogram). Note that PSD is
%the same as an autospectrum Sii.
%
% Inputs:
%
%       fs              [1 x 1]     Sampling rate in Hz
%       
%       data            [nS x nD]   Real (non-complex) sampled data in arbitrary
%                                   units. nD data series of length nS samples
%                                   can be done at the same time (independently)
%
% Outputs:
%
%       psd             [nF x nD]   Autospectrum (power spectral density) of the
%                                   data, one-sided at the nF frequencies in f
%
%       f               [nF x 1]    Frequencies in Hz at which the psd is
%                                   reported 
%
% Future Improvements:      none
%
% References:                                                                            none
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
% Revision History:        	11 December 2014    Created
%                           24 February 2015    Debugged to correct transposed
%                                               outputs, and issue of handling
%                                               only the first dimension of the
%                                               data. Minor update of header to
%                                               include copyright.
%                           10 March 2015       Got rid of the scaling by 1/fs
%
% Copyright (c) 2014-2015 Ocean Array Systems, All Rights Reserved.

% Number of samples. Quickest if this is a power of 2.
N = size(data,1);
if N/2 ~= round(N/2)
    % Zero pad by one element to make data an even length
    data = cat(1,data,zeros(1,size(data,2)));
    N = size(data,1);
end

% Take the magnitude of the FFT and scale it by fs*N
for i = 1:size(data,2)
    warning('COCKING AROUND - USE THE ONE WITHOUT fs AS PER MEYGEN DISCOVERY')
    psd(:,i) = (1/(fs*N)) * abs(fft(data(:,i))).^2; %#ok<AGROW>
%     psd(:,i) = (1/(N)) * abs(fft(data(:,i))).^2; %#ok<AGROW>
end

% Take the useful part of the FFT (real input; one sided)
psd = psd(1:N/2+1,:);

% Excluding DC and nyquist, scale by 2 to reflect the fact that you've taken one
% side only.
psd(2:end-1,:) = 2*psd(2:end-1,:);

% The frequencies for which PSD is reported are:
f = (0:fs/N:fs/2)';


end

