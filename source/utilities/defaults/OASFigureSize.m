function OASFigureSize(keyword, fh)
%OASFIGURES Function to set size of figure based on a keyword
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
% Copyright (c) 2014-2016 Ocean Array Systems, All Rights Reserved.

% Use current figure if no handle passed
if nargin == 1
    fh = gcf;
end

% Define figure sizes
switch keyword
    case 'onesquare'
        width  = 500;
        height = 450;
    case 'twosquare'
        width  = 1000;
        height = 450;
    case 'twostack'
        width = 500;
        height = 850;
        
    case 'threesquare'
        width  = 1500;
        height = 450;
    case 'onenarrow'
        width  = 500;
        height = 900;
    case 'twonarrow'
        width  = 1000;
        height = 900;
    case 'threenarrow'
        width  = 1500;
        height = 900;
    case 'onewide'
        width  = 1000;
        height = 450;
    case 'twowide'
        width  = 1000;
        height = 900;
    case 'threewide'
        width  = 1000;
        height = 1350;
    otherwise
        warning('Unknown keyword - define additional or set figure size manually')
end

% Ensure the top of the newly sized figure can't go off the screen by making the
% top right of the figure stay in the same place
current = get(fh, 'Position');
set(fh, 'Position', [current(1), (current(2)+current(4)-height), width, height])

