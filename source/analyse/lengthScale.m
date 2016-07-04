function [jLi, sigma, uvwBar, acorr] = lengthScale(dt, uvw, varargin)
%LENGTHSCALE Compute integral lengthscales of flow direct from data
%
% Syntax:
%       [jLi, sigma, uvwBar, acorr] = lengthScale(dt, uvw)
%           Uses eqs. 2,5 from reference [1] to directly calculate lengthscales
%           and standard deviations of the turbulent flow.
%
%       [jLi, sigma, uvwBar, acorr] = lengthScale(dt, uvw, plotTF)
%           Allows command line output and plot of autocorrelation functions to
%           be turned on or off (default on).
%
%       [jLi, sigma, uvwBar, acorr] = lengthScale(dt, uvw, plotTF, Tmax)
%           Allows command line output and plot of autocorrelation functions to
%           be turned on or off (default on).
%
% Inputs:
%
%       dt              [1 x 1]     Sampling period of u,v,w in seconds
%
%       uvw             [n x 3]     Speed in m/s of the flow sampled at a
%                                   constant rate (with period dt) where u is
%                                   streamwise, w is vertical upwards, and u,v,w
%                                   form a right handed set with v in the cross
%                                   stream direction.
%
%       plotTF          [1 x 1] boolean
%                                   True (default) causes a plot of the
%                                   autocorrelation peak to be made.
%
%       Tmax            [1 x 1]     Time in seconds beyond which there can be no
%                                   correlation. Typically improves accuracy by
%                                   removing white noise bias. Set to the
%                                   largest eddy turnover time.
%       
% Outputs:
%
%       jLi             [3 x 3]     [xLu xLv xLw; yLu yLv yLw; zLu zLv zLw]
%                                   Integral lengthscales of turbulence computed
%                                   by autocorrelation.
%                                   Warning; where mean flow speed in v,w
%                                   directions is zero, these integral
%                                   lengthscales are meaningless (and vanish)
%
%       sigma           [1 x 3]     Standard deviation of velocity fluctuations
%                                   in each of the three directions 
%                                       sigma = std(uvw-uvwBar)
%
%       uvwBar          [1 x 3]     Mean velocities in u,v,w directions produced
%                                   by the command >> mean(uvw,1) in m/s.
%
%       acorr           [2n-1 x 3]  Autocorrelation (normalised to unity via
%                                   xcorr's 'coeff' method) of each signal with
%                                   itself
%
% References:                                                                   
%
%   [1] M�cke T., Harkness C., and Argyriadis K., Offshore Wind Turbulence Model
%       vs. Measurements EWEA 2012 Presentation
%
%   [2] Garrad Hassan and Partners Document 282/BR/009 Issue 11 'GH Bladed
%       Theory Manual' July 2003
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
%                           24 February 2015    Updated inputs to include
%                                               plotting switch, allowing us to
%                                               turn off plots for batch running
%
% Copyright (c) 2014-2015 Ocean Array Systems, All Rights Reserved.

% Set to false to plot for debug
dbPlot = true;
if (nargin > 2) 
    if (islogical(varargin{1}))
        dbPlot = varargin{1}(1);
    else
        error('MATLAB:lengthscale:InvalidInput','Plotting switch must be [1 x 1] logical')
    end
end

% Filter out incoherent components to improve conditioning of the estimation
uvw = extractCoherent(uvw);

uvwBar = mean(uvw,1);
uvw = bsxfun(@minus,uvw,uvwBar);
% Autocorrelations
[acu] = xcorr(uvw(:,1),uvw(:,1),'coeff');
[acv] = xcorr(uvw(:,2),uvw(:,2),'coeff');
[acw] = xcorr(uvw(:,3),uvw(:,3),'coeff');
acorr = [acu(:) acv(:) acw(:)];

% Lose the first half (retaining the unity coefficient) and the last quarter 
startInd = ceil(numel(acu)/2);
endInd = ceil(3*numel(acu)/4);
acu = acu(startInd:endInd);
acv = acv(startInd:endInd);
acw = acw(startInd:endInd);

% Remove data from outside the viable range
if nargin > 3
    Tmax = varargin{2};
    nElements = floor(Tmax/dt);
    if nElements < numel(acu)
        acu = acu(1:nElements);
        acv = acv(1:nElements);
        acw = acw(1:nElements);
    end
end

% Integrate per ref [1] eq. 4 
Tuvw = [dt*trapz(acu) dt*trapz(acv) dt*trapz(acw)];

% Get the mean velocities in each direction
% uvwBar = mean(uvw,1);

% Ascertain the lengthscale 3x3 matrix by matrix multiplication
jLi = uvwBar'*Tuvw;

% Get the standard deviation of the flow
sigma = std(bsxfun(@minus, uvw, uvwBar));

if dbPlot
    disp('________________________________________________________________________')
    disp(' ')
    disp('Length and timescales, directly computed using autocorrelation')
    printScales(uvwBar, jLi, sigma, Tuvw)
    
    raiseFigure('lengthscale(): Autocorrelations acu, acv, acw');
    clf
    tau = (0:numel(acu)-1)*dt;
    plot(tau, acu,'b-')
    hold on
    plot(tau, acv,'g-')
    plot(tau, acw,'k-')
    xlabel('Delay time \tau (s)')
    ylabel('Normalised Autocorrelation')
    legend('uu','vv','ww')
    grid on
end

