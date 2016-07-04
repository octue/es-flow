function printScales(uvwBar, jLi, sigma, varargin)

disp('________________________________________________________________________')
disp(' ')
disp('Mean velocities (m/s) in x,y,z directions:')
disp(uvwBar)
disp(' ')
disp('Integral Lengthscales ordered as:    [xLu xLv xLw;')
disp('                                      yLu yLv yLw;')
disp('                                      zLu zLv zLw]')
disp(jLi)
disp(' ')
disp('Standard deviations sigma = std(uvw-uvwBar):')
disp(sigma)
disp(' ')
disp('Turbulent Intensities (assuming flow in u direction):')
disp(sigma/uvwBar(1))
disp(' ')
if nargin > 3
    disp('Integral Timescales (s) from u,v,w velocities:')
    Tuvw = varargin{1};
    disp(Tuvw)
end
disp('________________________________________________________________________')

