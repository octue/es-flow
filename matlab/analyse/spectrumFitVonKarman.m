function [jLiOpt, sigmaOpt, sigmaNoiseOpt] = spectrumFitVonKarman(freq, siiData, uvwBar, freqLimits, includeNoise, varargin)
%SPECTRUMFITVONKARMAN Nonlinear unbounded fit to determine parameters of 
% von Karman model that best represents a measured autospectrum.
%
% Syntax:
%       [jLiData, sigmaOpt, sigmaNoiseOpt] = spectrumFitVonKarman(dt, uvw)
%
%       [...] = spectrumFitVonKarman(dt, uvw, freqLimits)
%
%       [...] = spectrumFitVonKarman(dt, uvw, freqLimits, includeNoise)
%
%       [...] = spectrumFitVonKarman(dt, uvw, freqLimits, includeNoise, weights)
%
% Inputs:
%
%       f               [1 x 1]     Frequencies in Hz at which asd is given
%
%       asd             [n x 3]     Autospectral density S_uu, s_vv, S_ww
%
%       freqLimits      [1 x 2]     Extents of the frequency in Hz at which the
%                                   ASD is to be fitted. Must be real and
%                                   positive as [minf maxf]. A warning is issued
%                                   if maxf exceeds the nyquist criterion of the
%                                   data being fitted, in which case the fit
%                                   continues using the nyquist frequency as the
%                                   maximum limit.
%
%       includeNoise    [1 x 1]     Boolean. True causes standard deviation of
%                                   white noise to be included as a parameter in
%                                   the fit for each of the three directions.
%
%       weights         [n x 3]     Weightings in the range [0, 1] applied in
%                                   the least squares error estimator to the
%                                   different frequencies. This allows the user
%                                   to reduce the confidence in certain ranges
%                                   (e.g. if there is a greater error in the
%                                   data at higher frequencies, reduce the
%                                   weighting accordingly). The three columns
%                                   represent error in streamwise, cross stream
%                                   and vertical directions (x,y,z). As per the
%                                   behaviour of lscov(), Weights are 
%                                   applied outside the squared calculation so
%                                   weights*(differences^2) is minimised,
%                                   NOT (weights*differences)^2. This is
%                                   important where weights are calculated using
%                                   a physical model or reasoning.
%                                   Weights are applied for the nF frequencies
%                                   f that are output using the command
%                                       [~, f] = spectrumData(
%
% Outputs:
%
% References:
%
%   [1] RSI Internal Technical Note IN-037 'TiME, Islay Sound 2014
%       Nemo-MicroRider-Vector Measurements Synopsis', Rolf Lueck, 15 July 2014

%   [2] Garrad Hassan and Partners Document 282/BR/009 Issue 11 'GH Bladed
%       Theory Manual' July 2003
%
%   [3] Mücke T., Harkness C. and Argyriadis K. (2012) Offshore Wind Turbulence
%       Model vs. Measurement, EWEA Proceedings 2012
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
% Revision History:        	2 December 2014     Created
%                           8 December 2014     Documented and filled out the
%                                               function from the empty template



% At present we use manual weightings to determine the importance of different
% parts of the spectra in the different directions
if nargin > 5
    weights = varargin{1};
else
    % Uniformly weighted
    weights = ones(numel(freq),3);
end

% Apply frequency limits via weightings
weights(freq<min(freqLimits),:) = 0;
weights(freq>max(freqLimits),:) = 0;

% Define the first guess.
if includeNoise
    
    % First Guess
    %         jLiOpt = [5 4 3; 0 0 0; 0 0 0];
    %         TIOpt = [0.12 0.14 0.08];
    %         sigmaOpt = TIOpt*uvwBar(1);
    %         sigmaNoiseOpt = sqrt([0.005 0.0035 0.0007]);
    if nargin > 6
        x0 = varargin{2};
    else
        x0 = [5 4 3 0.2605 0.3040 0.1737 0.0707 0.0592 0.0265];
    end
    
else
    if nargin > 6
        x0 = varargin{2};
    else
        x0 = [5 4 3 0.2605 0.3040 0.1737];
    end
end

% Set up the per-iteration plots
raiseFigure('spectrumFitVonKarman() - Iteration Progress')
clf
subplot(2,2,1)
iter = 0;
lh1 = plot(iter,x0(1),'b-');
hold on
lh2 = plot(iter,x0(2),'g-');
lh3 = plot(iter,x0(3),'m-');
xlabel('Iteration #')
ylabel('Integral Lengthscale (m)')
grid on
legend('xLu','xLv','xLw')

subplot(2,2,2)
lh4 = plot(iter,x0(4),'b-');
hold on
lh5 = plot(iter,x0(5),'g-');
lh6 = plot(iter,x0(6),'m-');
xlabel('Iteration #')
ylabel('Standard Deviation \sigma_i')
grid on
legend('\sigma_u','\sigma_v','\sigma_w')

if includeNoise
    subplot(2,2,3)
    lh7 = plot(iter,x0(7),'b-');
    hold on
    lh8 = plot(iter,x0(8),'g-');
    lh9 = plot(iter,x0(9),'m-');
    xlabel('Iteration #')
    ylabel('Noise Standard Deviation \sigma_N_i')
    grid on
    legend('noise_u','noise_v','noise_w')
end

subplot(2,2,4)
lh10 = plot(iter,NaN,'b-');
hold on
lh11 = plot(iter,NaN,'g-');
lh12 = plot(iter,NaN,'m-');
lh13 = plot(iter,NaN,'r-');
xlabel('Iteration #')
ylabel('Least Squares Error between VK fit and modelled data')
grid on
legend('LSQ Error u','LSQ Error v','LSQ Error w','LSQ Error Sum')
    
    
% Set options for the search function
options = optimset(@fminsearch);
options = optimset(options,'MaxIter',10000,'MaxFunEvals',1e8);

% Perform the search
iterCtr = 0;
x = fminsearch(@nestEvalSpectrumFit, x0, options);
% if includeNoise
%     iterCtr = 0;
%     component = 1;
%     xu = fminsearch(@nestEvalSingleSpectrumFit, x0([1 4 7]), options);
%     iterCtr = 0;
%     component = 2;
%     xv = fminsearch(@nestEvalSingleSpectrumFit, x0([2 5 8]), options);
%     iterCtr = 0;
%     component = 3;
%     xw = fminsearch(@nestEvalSingleSpectrumFit, x0([3 6 9]), options);
% else
%     iterCtr = 0;
%     component = 1;
%     xu = fminsearch(@nestEvalSingleSpectrumFit, x0([1 4]), options);
%     iterCtr = 0;
%     component = 2;
%     xv = fminsearch(@nestEvalSingleSpectrumFit, x0([2 5]), options);
%     iterCtr = 0;
%     component = 3;
%     xw = fminsearch(@nestEvalSingleSpectrumFit, x0([3 6]), options);
% end

% Return results
% jLiOpt = [xu(1) xv(1) xw(1); 0 0 0; 0 0 0];
% sigmaOpt = [xu(2) xv(2) xw(2)]';
jLiOpt = [x(1:3); 0 0 0; 0 0 0];
sigmaOpt = x(4:6)';

if includeNoise
%     sigmaNoiseOpt = [xu(3) xv(3) xw(3)];
    sigmaNoiseOpt = x(7:9)';
else
    sigmaNoiseOpt = [0 0 0];
end


function sqDiff = nestEvalSpectrumFit(xParams)
    
    if any(xParams < 0)
        sqDiff = 10000;
        return
    end
    
    % Get the lengthscales from the input parameters
    % NB the vk model doesn't use yLi or zLi so zeros are fine but we just keep
    %    to the standard input form for this case
    jLi = [xParams(1:3); 0 0 0; 0 0 0];
    
    % All standard deviations to be absolute (negative values meaningless)
    xParams(4:end) = abs(xParams(4:end));
    
    % Get the standard deviations
    sigma = xParams(4:6)';
    
    % uvwBar and n already defined before the invocation of this nested function
    % and are in a shared workspace with the main function
    
    % Get the unnormalised autospectrum
    [siiVK] = spectrumVonKarmanBasic(jLi, uvwBar, freq, sigma);
    
    % If simulated noise is to be added as a parameter:
    %   From:
    %       http://www.gaussianwaves.com/2013/11/simulation-and-analysis-of-white-noise-in-matlab/
    %   [...] the power spectral density of the weakly defined white noise
    %   process is constant (flat) across the entire frequency spectrum. The
    %   value of the constant is equal to the variance or power of the white
    %   noise
    if numel(xParams) > 6
        [siiNoise] = repmat(xParams(7:9).^2, [numel(freq), 1]);
    else
        siiNoise = zeros(size(siiVK));
    end
    
    % siiData carries the unnormalised autospectrum directly computed from the
    % data for the frequencies given. We are optimising to minimise the least
    % squares difference between the two. We manually weight the importance of
    % the frequency components (see input help)
    sqDiffUVW = sum(weights.*(siiData-siiNoise-siiVK).^2,1);
    sqDiff = sum(sqDiffUVW);
    
    % Update the plots
    iterCtr = iterCtr + 1;
    set(lh1,'xData',0:iterCtr,'yData',[get(lh1,'yData'), xParams(1)])
    set(lh2,'xData',0:iterCtr,'yData',[get(lh2,'yData'), xParams(2)])
    set(lh3,'xData',0:iterCtr,'yData',[get(lh3,'yData'), xParams(3)])
    set(lh4,'xData',0:iterCtr,'yData',[get(lh4,'yData'), xParams(4)])
    set(lh5,'xData',0:iterCtr,'yData',[get(lh5,'yData'), xParams(5)])
    set(lh6,'xData',0:iterCtr,'yData',[get(lh6,'yData'), xParams(6)])
    if includeNoise
        set(lh7,'xData',0:iterCtr,'yData',[get(lh7,'yData'), xParams(7)])
        set(lh8,'xData',0:iterCtr,'yData',[get(lh8,'yData'), xParams(8)])
        set(lh9,'xData',0:iterCtr,'yData',[get(lh9,'yData'), xParams(9)])
    end
    set(lh10,'xData',0:iterCtr,'yData',[get(lh10,'yData'), sqDiffUVW(1)])
    set(lh11,'xData',0:iterCtr,'yData',[get(lh11,'yData'), sqDiffUVW(2)])
    set(lh12,'xData',0:iterCtr,'yData',[get(lh12,'yData'), sqDiffUVW(3)])
    set(lh13,'xData',0:iterCtr,'yData',[get(lh13,'yData'), sqDiff])
    
end % End nested function nestEvalSpectrumFit

function sqDiff = nestEvalSingleSpectrumFit(xParams)
    
    if any(xParams < 0)
        sqDiff = 10000;
        return
    end
    
    % Get the lengthscales from the input parameters
    % NB the vk model doesn't use yLi or zLi so zeros are fine but we just keep
    %    to the standard input form for this case
    switch component
        case 1
            
            % Lengthscales
            jLi = [xParams(1) 0 0; 0 0 0; 0 0 0];
    
            % Get the standard deviations
            sigma = [xParams(2) 0 0]';
            
            % Get the unnormalised autospectrum
            [siiVK] = spectrumVonKarmanBasic(jLi, uvwBar, freq, sigma);
            
        case 2
            
            % Lengthscales
            jLi = [0 xParams(1) 0; 0 0 0; 0 0 0];
    
            % Get the standard deviations
            sigma = [0 xParams(2) 0]';
            
            % Get the unnormalised autospectrum
            [siiVK] = spectrumVonKarmanBasic(jLi, uvwBar, freq, sigma);
            
        case 3
            
            % Lengthscales
            jLi = [0 0 xParams(1); 0 0 0; 0 0 0];
    
            % Get the standard deviations
            sigma = [0 0 xParams(2)]';
            
            % Get the unnormalised autospectrum
            [siiVK] = spectrumVonKarmanBasic(jLi, uvwBar, freq, sigma);
            
    end
    siiVK = siiVK(:,component);
    
    % If simulated noise is to be added as a parameter:
    %   From:
    %       http://www.gaussianwaves.com/2013/11/simulation-and-analysis-of-white-noise-in-matlab/
    %   [...] the power spectral density of the weakly defined white noise
    %   process is constant (flat) across the entire frequency spectrum. The
    %   value of the constant is equal to the variance or power of the white
    %   noise
    if numel(xParams) > 2
        siiNoise = (xParams(3).^2)*ones(size(siiVK));
    else
        siiNoise = zeros(size(siiVK));
    end
    
    % siiData carries the unnormalised autospectrum directly computed from the
    % data for the frequencies given. We are optimising to minimise the least
    % squares difference between the two. We manually weight the importance of
    % the frequency components (see input help)
    dfs = sqrt(abs((siiData(:,component)-siiNoise-siiVK)));
    sqDiff = sum(weights(:,component).*dfs);
    
    % Update the plots
    iterCtr = iterCtr + 1;
    switch component
        case 1
            set(lh1,'xData',0:iterCtr,'yData',[get(lh1,'yData'), xParams(1)])
            set(lh4,'xData',0:iterCtr,'yData',[get(lh4,'yData'), xParams(2)])
            if includeNoise
                set(lh7,'xData',0:iterCtr,'yData',[get(lh7,'yData'), xParams(3)])
            end
            set(lh10,'xData',0:iterCtr,'yData',[get(lh10,'yData'), sqDiff])
        case 2
            set(lh2,'xData',0:iterCtr,'yData',[get(lh2,'yData'), xParams(1)])
            set(lh5,'xData',0:iterCtr,'yData',[get(lh5,'yData'), xParams(2)])
            if includeNoise
                set(lh8,'xData',0:iterCtr,'yData',[get(lh8,'yData'), xParams(3)])
            end
            set(lh11,'xData',0:iterCtr,'yData',[get(lh11,'yData'), sqDiff])
        case 3
            set(lh3,'xData',0:iterCtr,'yData',[get(lh3,'yData'), xParams(1)])
            set(lh6,'xData',0:iterCtr,'yData',[get(lh6,'yData'), xParams(2)])
            if includeNoise
                set(lh9,'xData',0:iterCtr,'yData',[get(lh9,'yData'), xParams(3)])
            end
            set(lh12,'xData',0:iterCtr,'yData',[get(lh12,'yData'), sqDiff])
    end
    
end % End nested function nestEvalSingleSpectrumFit

end % End main function

