function sigOut = wsola(sigIn, fs, scaleFactor, varargin)
%WSOLA Rescales the time axis of a time series dataset
% Uses Waveform Similarity Overlap-Add (WSOLA) algorithm described in Verhelst
% and Roelands at:
%   http://www.etro.vub.ac.be/research/dssp/PUB_FILES/int_conf/ICASSP-1993.pdf
%
% The aim of a WSOLA time scaling algorithm is to produce a synthetic waveform
% y(n) that maintains maximal local similarity to the original waveform x(m) in
% corresponding neighbourhoods of related sample indices n = tau(m).
%
% Local similarity is defined somewhat vaguely using a 'similarity' metric,
% which can be based on cross correlation or cross 
%
% Syntax:  
%
%       sigOut = wsola(sigIn, fs, scale) Rescales a time series signal to a new
%       length at the same sampling frequency, using the WSOLA algorithm.
%
%       sigOut = wsola(sigIn, fs, scale, method) allows selection of the method
%       used to provide a similarity metric. Available methods are 'xcorr'
%       (normalised cross correlation) or the default 'amdf' (cross Average
%       Magnitude Difference Function).
%
%       sigOut = wsola(sigIn, fs, scale, method, windowSize) selects window size
%       in terms of number of samples. Default is 200, which corresponds to 20ms
%       window size for a 10kHz signal.
%
%       sigOut = wsola(sigIn, fs, scale, method, windowSize, overlap) selects
%       the overlap ratio where 0 <= overlap < 1 and the default is 0.5 (a 50%
%       window overlap).
%
% Inputs:
%    
%       sigIn               [1 x ns]    Data signal sampled monotonically in
%                                       time at frequency fs for a period ns/fs
%
%       fs                  [1 x 1]     Sampling frequency in Hz
%
%       scale               [1 x 1]     Scale factor by which timeseries must be
%                                       lengthened.
%
%       method              string      Specifies the method used to provide
%                                       a similarity measure for the WSOLA
%                                       algorithm, either 'xcorr' (normalised
%                                       cross correlation) or the default 'amdf'
%                                       (cross average magnitude difference
%                                       function). See figure (3) of the
%                                       reference for more information.
%
%       windowSize          [1 x 1]     The window size expressed as an integer
%                                       number of samples. To express windowsize
%                                       in terms of time, use:
%                                           windowSize = windowTime*fs
%
%       overlap             [1 x 1]     0<= overlap < 1
%                                       The overlap of successive windows used
%                                       for the rescaling process. A value of
%                                       0.5 (50%) is the most common choice,
%                                       balancing signal continuity and speed.
%
% Outputs:                  
%
%       sigOut              [1 x round(ns*scale*)]
%                                       Data signal sampled monotonically in
%                                       time at a frequency fs, stretched over
%                                       time by the input scale factor without
%                                       modifying the signal pitch.
%
% References:
%       [1] www.etro.vub.ac.be/research/dssp/PUB_FILES/int_conf/ICASSP-1993.pdf
%
% Future Improvements:      none
%
% Author:                   T. H. Clark
% Work address:             3 Charles Babbage Road
%                           Cambridge
%                           CB3 0GT
% Email:                    tom.clark@oceanarraysystems.com
% Website:                  www.oceanarraysystems.com
%
% Revision History:         13 August 2014      Created

% Main body of code and subfunctions copyright:
% Copyright (c) 2013, Arthur Dgn
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% Defaults
similMethod = 'amdf';
windowSize = 200;
overlap = 0.5;

% Parse nondefault inputs and check validity
if nargin > 3
	similMethod = lower(varargin{1});
    if (~ischar(similMethod)) || (~any(strcmpi(similMethod,{'xcorr','amdf'})))
        error('MATLAB:InvalidInput','wsola: Input ''method'' must be a string describing the similarity method, either ''xcorr'' or ''amdf'' are valid')
    end
end
if nargin > 4
    windowSize = varargin{2};
    if (~isnumeric(windowSize)) || (windowSize ~= round(windowSize)) || (windowSize < 1)
        error('MATLAB:InvalidInput','wsola(): Input ''windowSize'' must be a positive integer value')
    end
end
if nargin > 5
    overlap = varargin{3};
    if (~isnumeric(overlap)) || (overlap < 0) || (overlap >= 1)
        error('MATLAB:InvalidInput','wsola(): Input ''overlap'' must be a numeric value in the range 0 <= overlap < 1')
    end
end

error('TOM BE VERY CAREFUL - THIS NEVER WORKED PROPERLY (check the end of the signal)')

sf_split_margin = 3;
orig_len = length(sigIn);
% recursively treat scale factors too big or too small
if (scaleFactor > sf_split_margin) || (scaleFactor < 1/sf_split_margin)
    num_recur = ceil(abs((log(scaleFactor)./log(sf_split_margin))));
    
    scaleFactor_recur = scaleFactor.^((num_recur-1)/num_recur);
    
    sigIn = wsola(sigIn, fs, scaleFactor_recur, varargin{:});
    
    scaleFactor = scaleFactor./(length(sigIn)./orig_len);
end


% Get constants used in original code from the parsed and checked inputs
win_time = windowSize/fs;
overlap_ratio = overlap;
max_err = min(0.005, win_time*overlap_ratio./scaleFactor/2);

if max_err == 0.005
    warning('AUTHOR NOTE - TOM - NOTE YOU ARE ON THE MAX ERR LIMIT. FIND OUT WHAT IT DOES')
end

% lengths
win_len = ceil(win_time.*fs);
max_err_len = ceil(max_err.*fs);
step_len = floor(overlap_ratio.*win_len);

% vectors
win = hann(win_len);
% orig_scale = 1:length(sigIn);
new_scale = (1:round((length(sigIn).*scaleFactor*2)))';
sigOut = zeros(size(new_scale));
win_out = zeros(size(new_scale));

cursor_in = 1;
cursor_out = 1;

while cursor_in<(length(sigIn)-win_len-max(step_len, step_len./scaleFactor) - 2*max_err_len)...
                                    && cursor_out<(length(sigOut)-win_len)
    % input segments
        size(sigIn(cursor_in:(cursor_in+win_len-1)))
        size(win)
    new_seg = sigIn(cursor_in:(cursor_in+win_len-1)).*win;
    new_seg_neihbour = sigIn((cursor_in+step_len):(cursor_in+step_len+win_len-1)).*win;
    % overlap add
    sigOut(cursor_out:(cursor_out+win_len-1)) = ...
            sigOut(cursor_out:(cursor_out+win_len-1))+...
            new_seg;
    % overlap add window normalization vec
    win_out(cursor_out:(cursor_out+win_len-1)) = ...
            win_out(cursor_out:(cursor_out+win_len-1))+...
            win;
        
    % move cursors
    cursor_out = cursor_out + step_len;
    cursor_in = cursor_in + round(step_len./scaleFactor);
    % new candidate
    new_seg_cand = sigIn((cursor_in):(cursor_in+win_len-1)).*win;
    % similarity calc. Note: matlab apears to be using FFT for xcorr, and
    % so computes the whole thing (instead of just the center). 
    % Writing own version of xcorr would probably make it faster.
    if strcmp(similMethod, 'xcorr')
        shift = max_xcorr_similarity(new_seg_neihbour, new_seg_cand, max_err_len);
    elseif strcmp(similMethod, 'amdf')
        shift = min_amdf_similarity(new_seg_neihbour, new_seg_cand, max_err_len);
    else
        return
    end 
    % adjust cursor place
    cursor_in = cursor_in - shift;
end

% remove slack
sigOut(cursor_out:end) = [];
win_out(cursor_out:end) = [];

sigOut = sigOut./(win_out+eps); %normalize to remove possible modulations

end

    % cross correlation based similarity calculation
    function shift = max_xcorr_similarity(seg1, seg2, max_lag)
        [~, max_i] = max(xcorr(seg1, seg2, max_lag, 'unbiased'));
        shift = max_i - max_lag;
    end

    % amdf based similarity calculation
    function shift = min_amdf_similarity(seg1, seg2, max_lag)
        n = length(seg1);
        amdf = ones(1,2*max_lag-1);
        for lag=-max_lag:max_lag
            amdf(lag+max_lag+1) = sum(abs(seg2(max(1,(-lag+1)):min(n,(n-lag)))-...
                           seg1(max(1,(lag+1)):min(n,(n+lag))) ))/n;
        end
        [~, min_i] = min(amdf);
        shift = min_i - max_lag;
    end