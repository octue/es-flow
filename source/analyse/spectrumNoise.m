function [siiNoise] = spectrumNoise(freq, Nsamp, fs, meanValue, stdDev)
%SPECTRUMNOISE Outputs the autospectrum of white (gaussian) noise applied around
%a mean value
%
% Syntax:
%
% Inputs:
%
%       freq            [nf x 1]    Frequency in Hz at which the autospectrum is
%                                   output.
%
%       Nsamp           [1 x 1]     The number of data points used to create the
%                                   artificial spectrum (should be the same as
%                                   the number of data points in the spectrum
%                                   you're comparing it to or removing it from).
%
%       fs              [1 x 1]     Sampling frequency (Hz) used to create the
%                                   artificial autospectrum (should be the same
%                                   as the number of data points in the spectrum
%                                   you're comparing it to or removing it from).
%
%       meanValue       [1 x nd]    Mean value around which gaussian noise is
%                                   applied.
%
%       stdDev          [1 x nd]    Standard deviation of the gaussian noise.
%
% Outputs:
%
%       siiNoise        [nf x nd]   Autospectra of the simulated noisy signals
%
% Note:
%
%   [1] Considered using two values of stdDev for each signal; where one gives
%       an absolute error and one is scaled proportional to the mean value.
%       Whilst somewhat useful if simulating noise from an instrument this is
%       easily achieved by addition of the standard deviations required on the
%       input; the two sources of error, provided they are both white, have the
%       same effect on the shape of the spectrum and therefore cannot be
%       separated / identified independently from the characteristics of a noisy
%       signal.
%
% References:
%
%   [1] 
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
% Revision History:        	8 December 2014     Created
%                           9 December 2014     Implemented
%                           10 December 2014    Got FUCKED OFF with the FUCKING
%                                               DC component scaling and totally
%                                               reimplemented as analytical
%                                               determination of the spectrum

% Input Checks
% 1. Check if meanValue and stdDev are both 1 x nDims
% 2. Check stdDev > 0
% 3. Check freq is real and positive and correct shape (column vector)
% 4. Check N and fs are scalar

% Compute non-dc power levels Using Jean Baptiste Richard's EWTEC paper
siiNoise = repmat(2*stdDev.^2./fs,[numel(freq),1]);

% Compute power at the DC level and place it into any zero frequency entries.
% Clumsy code using repmat...
dc = repmat((Nsamp/fs)*meanValue.^2,[numel(freq),1]);
siiNoise(freq==0,:) = dc(freq==0,:);

