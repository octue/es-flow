function [varargout] = raiseFigure(figTitle)
%RAISEFIGURE Creates/raises a figure by title and tag, rather than handle/number
% Replaces MATLAB's
%           fh = figure()   and
%                figure(h)
% commands where a user wishes to identify a figure by its title (and tag), 
% rather than by its handle.
%
% This is useful in cases where the same figure needs to be accessed by multiple
% workspaces.
%
% Example:
%   A function (which is called multiple times) updates a figure each time
%   (displaying some data). There are many ways to do this, but typically, the 
%   plot must be initialised before the first function call, and the figure
%   handle passed with each call to the function.
%
%   Instead, raiseFigure('FigTitle') can be used within the function call
%   (instead of figure() ). Each time that function is executed, it'll either
%   create or make current the same figure (which has the title 'FigTitle'):
%
%   function [] = exampleFunction(args)
%
%       fh = raiseFigure('My Test Figure')
%       
%       % Clear it (if you want to)
%       clf
%
%       % Set some properties and/or plot some data
%       set(fh,'someProperty',someValue)
%       plot(someData)
%
%
%
% Background:
%   I coded this because I had functions with figure numbers hard-coded. So I'd
%   call figure(801) (or whatever number) and each time the function was called,
%   it'd plot into figure 801. Of course, sooner or later I used the same number
%   twice and got into a massive mess!
%
% Syntax:
%            raiseFigure(figTitle)
%               Creates or makes current a figure whose tag corresponds to the
%               string figTitle. The figure title is also set to
%               that string, although may be subsequently changed.
%
%       fh = raiseFigure(figTitle)
%               As per raiseFigure(figTitle). Returns the numeric handle to
%               the figure
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
% Copyright (c) 2014 Ocean Array Systems, All Rights Reserved.

% Check that the input argument is a character string
if (~ischar(figTitle)) || (size(figTitle,1) ~= 1)
    error('MATLAB:raiseFigure:notAString','Figure title (tag) input must be a character string, e.g. ''MyTitle''')
end
    
% Create (or raise existing) figure
found_fh = findobj('Tag',figTitle);
if isempty(found_fh)
    % Create
    fh = figure('Name',figTitle,'Tag',figTitle, 'NumberTitle','off');
else
    % Make the figure current (and return handle)
    fh = figure(found_fh(1));
end

% Issue a warning if more than one figure exists with that tag
if numel(found_fh) > 1
    warning('MATLAB:raiseFigure:multipleFiguresDetected','More than one figure with the same tag was detected. The incorrect figure may have been raised')
end

% If called with arguments
if nargout > 0
    varargout{1} = fh;
end
    
    


