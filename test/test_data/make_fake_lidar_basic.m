function make_fake_lidar_data(file)
%MAKE_FAKE_LIDAR_DATA Saves a fake LiDAR dataset to a .mat file
%
% Syntax:  
%       make_fake_lidar_data(filename)
%
% Inputs:
%
%       file        string      filename to save to. Will be overwritten if it
%                               exists already. The path must already exist.
%
% References:
%
%   [1] Perry AE and Marusic I (1995) A wall-wake model for turbulent boundary
%       layers. Part 1. Extension of the attached eddy hypothesis J Fluid Mech
%       vol 298 pp 361-388
%
%   [2] I. Marusic and A. E. Perry, ?A wall-wake model for the turbulence 
%       structure of boundary layers. Part 2. Further experimental support.,? 
%       J. Fluid Mech., vol. 298, 1995.
%
%   [3] H. Oertel, Ed., Prandtl?s Essentials of Fluid Mechanics, Second. Ed. 
%       2003.
%
% Future Improvements:  
%
%   [1] Add deliberate corruptions to the dataset (time discontinuities, NaNs
%   etc).
%
%   [2] Add confidence metrics per real LiDAR values.
%
%   [3] Add Ekman effect with realistic variation of turbulent diffusion to
%       to give realistic deviation (~20degrees) at ground level
%
% Author:                   T. Clark
% Work address:             Ocean Array Systems Ltd
%                           Hauser Forum
%                           3 Charles Babbage Road
%                           Cambridge
%                           CB3 0GT
% Email:                    tom.clark@oceanarraysystems.com
% Website:                  www.oceanarraysystems.com
%
% Copyright (c) 2016-17 Ocean Array Systems, All Rights Reserved.

% LiDAR, sampling 100 bins at 1m spacing, with half angle of 5 degrees, for 24
% hours at 8 Hz, with the first measurement volume at 5m
type = 'lidar_basic';
n_samples = 24*60*60*8;
n_bins = 100;
half_angle = 5;
z = ((0:1:99) + 5)';
u = (rand(n_bins, n_samples)-0.5)*0.1;
v = (rand(n_bins, n_samples)-0.5)*0.1;
w = (rand(n_bins, n_samples)-0.5)*0.1;
t_seconds = linspace(0, (n_samples-1)*0.125, n_samples)';
t_start_vec = [2014 12 10 13 41 23.2540];
t = 86400*(datenum(t_start_vec) - datenum([1970 01 01 0 0 0])) + t_seconds;
position = [52.227799, 0.263672];

units.half_angle = 'degrees';
units.z = 'm';
units.position = 'wgs84_degrees';
units.u = 'm/s';
units.v = 'm/s';
units.w = 'm/s';
units.t = 'posix_s';

% Assume a mean velocity profile based on theory of Ref [1] and values of Pi and
% S from Ref [2], the fully developed case.
Pi = 3.23;
S = 38.4;

% Assume atmospheric boundary layer thickness of order 1km 
deltac = 1000;
U1 = 20;
[U_mag] = getMeanProfile(Pi, S, deltac, U1, z);

% Assume turbulent visosity for the atmosphere as per Ref [3] section 12.1.5
% (p.582)
% nu_tau = 10; % NOT USED YET - SEE ENHANCEMENT [3]

% Apply an Ekman variation like 1/z, normalised to give us a generic shape with
% a 20 degree deviation at 5m height. Could be much more realistic.
inv_ht = 1./z;
inv_ht_norm = (inv_ht - 1/deltac)./(max(inv_ht) - 1/deltac);
ekman_angle = 20*inv_ht_norm;
U_x = U_mag.*cosd(ekman_angle);
U_y = U_mag.*sind(ekman_angle);

% Apply the Ekman varied mean flow profile to the base fluctuations and scale
% the fluctuations by U1
u = bsxfun(@plus, U_x, u*U1);
v = bsxfun(@plus, U_y, v*U1);
w = w*U1;

% % Plot a check figure
% raiseFigure('Quiver check')
% clf
% xy = zeros(size(z));
% for i = 1:20;
%     quiver3(xy, xy, z, u(:,i), v(:,i), w(:,i))
%     hold on
% end

% Save the test file
save(file, 't', 'z', 'u', 'v', 'w', 'position', 'half_angle', 'type', 'units', '-v7.3')




