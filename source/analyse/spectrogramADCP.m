function [S F T Pxx] = spectrogramADCP(adcpData, bin, U, varargin)
%SPECTROGRAMADCP Time dependent spectrogram of velocity content at a bin
%   Allows a quick and easy view of where in a tidal cycle the frequency
%   content alters substantially. This is useful as a first look in turbulence
%   analysis, as it tells us where within the tidal cycle there are regime
%   changes that fundamentally affect the turbulent content of the flow.
%
% Syntax:  
%
%       [S F T] = spectrogramADCP(adcpData, binIndex, U)
%                  Computes a spectrogram of the u velocity signal at
%                  bin(binIndex) in adcpData. Uses a hanning window with zero
%                  overlap set (usingthe adcpData time vector) to span 10
%                  minutes of fluid flow. 
%
%       [S F T] = spectrumADCP(..., 'windowMinutes', windowMinutes)
%                  Computes spectrum using a hanning window spanning
%                  windowMinutes minutes of the default 10 minutes.
%
%         [...] = spectrumADCP(..., 'Plot', true)
%                  Adds output to a figure.
%
% Inputs:
%
%       adcpData        structure       OAS standard adcp data structure,
%                                       containing time series velocity data.
%                                       See help('loadADCP') for fields
%                                       reference.
%
%       bin             [1 x 1]         Index of the bin of ADCP data from which
%                                       to take the PSD.
%
%       U               [1 x 1]         Reference velocity used for scaling
%                                       (in m/s)
%
% Optional Inputs:
%
%       Parameter       Value
%
%       windowMinutes   [1 x 1]         Time in minutes spanned by the window
%                                       over which each short time spectrum is
%                                       taken. Columns in output S correspond to
%                                       PSDs taken windowMinutes apart.
%
%       plot            [1 x 1]         Boolean, default false
%                                       Creates a figure detailing the
%                                       PSD with relevant overlays
%
%       overlap         [1 x 1]         Default 0 (recommended)
%                                       Overlap as a fraction of window size.
%                                       Must be 1 > overlap >= 0.
%           
%
% Outputs:
%
%       See help(spectrogram)
%
% References:
%
%       [1] Richard J.B., Thomson J., Polagye B., Bard J., (2013) Method for
%           identification of Doppler Noise Levels in Turbulent Flow
%           Measurements Dedicated to Tidal Energy. Proceesdings of EWTEC,
%           Aalborg, September 2013
%
% Future Improvements:      
%
%       [1] Addition of GPU based FFT functionality
%
% Other m-files required:   none
% Subfunctions:             none
% Nested functions:         none
% MAT-files required:       none
%
% Author:           T. H. Clark
% Work address:     Hauser Forum
%                   3 Charles Babbage Road
%                   Cambridge
%                   CB3 0GT
% Email:            tom.clark@oceanarraysystems.com
% Website:          www.oceanarraysystems.com
%
% Created:          23 July 2014
% Revisions:        

% Set default parameters
opts.windowMinutes = 10;
opts.plot   = false;
opts.useGPU = false;
opts.overlap = 0;
opts = parse_pv_pairs(opts,varargin);



%% CHECK ON INPUTS

% Check the time signal is monotonic

% Check overlap is 0 < overlap < 1. Error if >= 1 and warn of lost data if < 0



%% SIGNAL PROCESSING

% Get signal 1
u = adcpData.u(bin,:);
v = adcpData.v(bin,:);
w = adcpData.w(bin,:);

% Determine the overlap in terms of window length 
dt = adcpData.t(2) - adcpData.t(1);
fs = 1/dt;
windowSize = nextpow2(datenum([0 0 0 0 opts.windowMinutes 0])/dt);
windowDuration = datevec(dt*windowSize);
disp(['spectrogramADCP: Using nearest (next power of two elements) window duration of ' num2str(windowDuration(5)) ' minutes.']) 
noverlap = floor(windowSize*opts.overlap);
disp(['spectrogramADCP: Using overlap of ' num2str(noverlap)]) 

% Build the spectrogram
if opts.plot || (nargout > 3)
    [S, F, T, Pxx] = spectrogram(u, hamming(windowSize), noverlap, windowSize, fs);
else
    [S, F, T, Pxx] = spectrogram(u, hamming(windowSize), noverlap, windowSize, fs);
end



%% PLOT OUTPUTS

if opts.plot
    
    raiseFigure(['Spectrogram at height z = ' num2str(adcpData.z(bin))])
    clf
    
    % Plot the spectrogram with time along the x axis
    if length(T)==1
        % surf requires a matrix for the third input.
        args = {[0 F],F,10*log10(abs([Pxx Pxx])+eps)};
    else
        args = {T,F,10*log10(abs(Pxx)+eps)};
    end
    surf(args{:},'EdgeColor','none');
    axis xy; axis tight;
    view(0,90);
    ylabel('Frequency [Hz]')
    
    % Plot horizontal line denoting the beam separation frequency
    kbs = fbs*2*pi/U;
    plot([adcpData.t(1) adcpData.t(end)], [fbs fbs], 'r-')
    
    % Plot the input signal on axes directly below
    subplot(2,1,2)
    plot(adcpData.t,u,'b-')
    ylim(max(abs(u))*[-1.1 1.1])
    ylabel('Streamwise velocity u')
    datetick
    
    % Label the axes
    xlabel('Time of data acquisition')
    
    % Add legend
    legend('PSD(u)','Beam Separation Limit')
    
    % Add a colorbar for PSD scale
    colorbar
    
    
    
end













