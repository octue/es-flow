function [asdout] = spectrumVonKarmanBasic(jLi, uvwBar, f, varargin)
%SPECTRUMVONKARMANBASIC Outputs the Autospectral Density of the basic von Karman
%model in the frequency range specified.
%
% The Basic von Karman model is given by Reference [1] as:
%
% $$ \frac{n S_{uu}(n)}{{\sigma_u}^2} = \frac{ 4\tilde{n_u}} {(1+70.8\tilde{n_u}^2)^{5/6} } $$
%
% and
%
% $$ \frac{n S_{ii}(n)}{{\sigma_i}^2} = \frac{ 4\tilde{n_i}(1 + 755.2\tilde{n_i}^2)} {(1+282.3\tilde{n_i}^2)^{11/6} } $$
%
% where
% 
% $$ \tilde{n_i} = \frac{n ^xL_i}{\overbar{u}} \qquad \forall \qquad i = v,w $$
%
% Syntax:
%
%      [nasd] = spectrumVonKarmanBasic(jLi, uvwBar, f)
%               Computes the normalised autospectral density nS_i(n)/sigma_i^2
%               at frequencies f from turbulent lengthscales jLi and mean flow
%               velocities uBar, using the basic von Karman model.
%
%       [asd] = spectrumVonKarmanBasic(jLi, uvwBar, f, sigma)
%               Computes the un-normalised autospectral density S_u(n).
%
% Inputs:
%
%       jLi             [3 x 3]     [xLu xLv xLw; yLu yLv yLw; zLu zLv zLw]
%                                   Integral lengthscales of turbulence computed
%                                   by autocorrelation.
%                                   Warning; where mean flow speed in v,w
%                                   directions is zero, these integral
%                                   lengthscales are meaningless (and vanish)
%
%       uvwBar          [1 x 3]     Mean velocities of flow in u,v,w directions
%                                   in m/s.
%
%       f               [n x 1]     Frequencies (Hz) at which asd is calculated
%
%       sigma           [1 x 3]     Standard deviation of velocity fluctuations
%                                   in each of the three directions 
%                                       sigma = std(uvw-uvwBar)
%       
% Outputs:
%
%       asd             [n x 3]     Autospectral density S_ii for i = u,v,w
%
%       nasd            [n x 3]     Normalised autospectral density
%                                   f*S_ii/sigma_i^2
%
% References:
%
%   [1] Garrad Hassan and Partners Document 282/BR/009 Issue 11 'GH Bladed
%       Theory Manual' July 2003
%
%   [2] Veers, P. S. Three Dimensional Wind Simulation, SAND88 - 0152, Sandia
%       National Laboratories, March 1988
%
%   [3] Mücke T., Harkness C. and Argyriadis K. (2012) Offshore Wind Turbulence
%       Model vs. Measurement, EWEA Proceedings 2012
%
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
% Revision History:        	4 December 2014     Created

% Get values out of jLi in more familiar notation:
xLu = jLi(1,1);
xLv = jLi(1,2);
xLw = jLi(1,3);

% Get streamwise average velocity
u10 = uvwBar(1);

% Get the normalised frequencies
f = f(:); % Ensure columnar
tildeNu = (f.^xLu)/u10;
tildeNv = (f.^xLv)/u10;
tildeNw = (f.^xLw)/u10;

% Get the normalised spectrum
nasdu = 4*tildeNu                         ./ (1 +  70.8*tildeNu.^2).^( 5/6);
nasdv = 4*tildeNv.*(1 + 755.2*tildeNv.^2) ./ (1 + 282.3*tildeNv.^2).^(11/6);
nasdw = 4*tildeNw.*(1 + 755.2*tildeNw.^2) ./ (1 + 282.3*tildeNw.^2).^(11/6);


% If sigma is also given...
if nargin > 3
    
    % Denormalise the output autospectra
    sigma = varargin{1};
    asdu = nasdu*(sigma(1)^2)./f;
    asdv = nasdv*(sigma(2)^2)./f;
    asdw = nasdw*(sigma(3)^2)./f;
    
    % Fix the NaN bug at zero frequency
    asdu(isnan(asdu)) = 0;
    asdv(isnan(asdv)) = 0;
    asdw(isnan(asdw)) = 0;
    
    asdout = [asdu asdv asdw];
    
else
    
    % Output the normalised autospectra
    asdout = [nasdu nasdv nasdw];

end

