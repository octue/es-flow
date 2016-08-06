function [k1, psd, fbs, sbs, N, K, fcut] = spectrumADCP(adcpData, bin, U, varargin)
%SPECTRUMADCP Power Spectral Density (PSD) of velocity content at a bin
%
% Syntax:  
%
%       [k1 psd] = spectrumADCP(adcpData, binIndex, U)
%                  Determines the power spectral density of the velocity
%                  fluctuations within bin binIndex of adcpData.
%
%       [k1 psd fbs sbs] = spectrumADCP(adcpData, binIndex, U)
%                  Determines PSD and frequency of beam separation (defined in
%                  ref [1]) above which measurements of features in ADCP data
%                  are inherently corrupted. The corresponding scale is also
%                  returned, below which features are corrupted.
%
%       [k1 psd fbs sbs N K fcut] = spectrumADCP(adcpData, binIndex, U)
%                  Determines coefficients N, K and fcut relating to Doppler
%                  Noise Level as described in Reference [1].
%
%       [k1 psd fbs sbs N K fcut] = spectrumADCP(adcpData, binIndex, U, freqRange)
%                  Determines coefficients N, K and fcut relating to Doppler
%                  Noise Level as described in Reference [1] using only data
%                  within a specified frequency range (to be used with care).
%
%       [...] = spectrumADCP(..., 'Parameter', value)
%                  Provides custom parameters as listed below.
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
%       plot            [1 x 1]         Boolean, default false
%                                       Creates a figure detailing the
%                                       PSD with relevant overlays
%
% Outputs:
%
%       k1              [nf x 1]        Streamwise wavenumbers k1 = 2*pi*f/U
%
%       psd             [nf x 1]        Power spectral density of the flow based
%                                       on streamwise velocity component u.
%
%       fbs             [1 x 1]         Frequency of beam separation, above
%                                       which the signal is corrupted by the
%                                       inherent limitations of a beam-separated
%                                       ADCP measurement. This is defined in [1]
%                                       as fbs = U/z which relies on an
%                                       approximation that z ~= b.
%                                       Here, beam angle is taken from the adcp
%                                       data structure and fbs computed directly
%                                       as:
%                                       fbs = U/b where beam separation b = b(z,
%                                       theta) for beam separation angle theta.
%
%       sbs             [1 x 1]         The physical lengthscale in m
%                                       corresponding to fbs below which
%                                       physical features are unlikely to be
%                                       captured correctly by the ADCP.
%                                       sbs = b ( ~= z for beam angles 20-25
%                                       degrees).
%
%       N               [1 x 1]         Noise level coefficient from Eq. 10,
%                                       Reference [1]
%
%       K               [1 x 1]         PSD signal coefficient K from Eq. 10,
%                                       Reference [1]
%
%       fcut            [1 x 1]         Cutoff frequency above which noise
%                                       dominates the signal (see Eq. 11
%                                       Reference [1] for strict definition).
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
%       [2] Ability to weight the frequency range rather than simply provide a
%           range of valid frequencies for estimation of the doppler noise level
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
opts.plot   = false;
opts.useGPU = false;
if nargin > 3
    opts = parse_pv_pairs(opts,varargin(2:end));
end
if nargin > 2
    range = varargin{1};
else
    range = [-Inf Inf];
end


%% CHECK ON INPUTS

% Check the time signal is monotonic


%% SIGNAL PROCESSING

% Get signal 1
u = adcpData.u(bin,:);
v = adcpData.v(bin,:);
w = adcpData.w(bin,:);

% Base the FFT on flow speed, not u1
speed = sqrt(u.^2);
Y = fft(speed);

% Scale and magnitude of one-sided FFT to get PSD and frequencies
n=length(Y);
psd = abs(Y(1:floor(n/2))).^2;
dt = (adcpData.t(2) - adcpData.t(1));
nyquist = 1/(2*dt);
freq = (1:n/2)/(n/2)*nyquist;
k1 = 2*pi*freq.*U;

warning('Developer note: Check spectrumADCP() against the spectrumData() which was validated during a commercial work package, and if necessary rewrite to utilise spectrumData')

% Cheap computation of fbs and sbs
sbs = 2*tand(adcpData.beamAngle)*(adcpData.z(bin) - adcpData.zUnit);
fbs = U/sbs;



%% PERFORM DOPPLER NOISE LEVEL ESTIMATION

if nargout > 4
    
    % fit to a -5/3 gradient
    freq53 = freq.^(-5/3);
    
    % Weight to exclude undesired frequencies.
    ci = ones(size(freq53));
    ci(freq>max(range)) = 0;
    ci(freq<min(range)) = 0;
    
    % Assemble the matrices in Eq. 10 Reference [1]
    A = [sum(ci) sum(ci.*freq53); sum(ci.*freq53) sum(ci.*freq53.*freq53)];
    b = [sum(ci.*psd); sum(ci.*psd.*freq53)];
    
    % Solve
    x = A\b;
    N = x(1);
    K = x(2);
    
    % Get cutoff freq from Eq. 11 Reference [1]
    fcut = (N/K)^(-3/5);
    
end

    
    
%% PLOT OUTPUTS

if opts.plot
    
    raiseFigure(['Power Spectral Density at height z = ' num2str(adcpData.z(bin))])
    clf
    
    % Plot the PSD
    plot(k1,psd,'b-')
    hold on 
    
    % Plot vertical line denoting the beam separation frequency
    kbs = fbs*2*pi/U;
    plot([kbs kbs], [0 max(psd)],'r-')
    
    % Plot the lines denoting N and K
    if nargout > 4
        
        % Plot the vertical line denoting the cutoff frequency
        kcut = fcut*2*pi/U;
        plot([kcut kcut],[0 max(psd)],'r--')
        
        % Plot the gradient and intercept lines from the model fit
        plot(k1, K*freq53, 'm-')
        plot([min(k1) max(k1)], [N N], 'm--')
        
    end
        
    % Label the axes
    ylabel('PSD [(ms^-^1)^2.Hz^-^1]')
    xlabel('Streamwise wavenumber k1 = 2\pi f/\bar(U(z))')
    
    % Add legend
    if nargout > 4
        legend('PSD (flow speed)', 'Beam Separation Limit', 'Noise Cutoff Frequency', 'K*f^-^5^/^3', 'N')
    else
        legend('PSD (flow speed)','Beam Separation Limit')
    end
    
    
    
end













