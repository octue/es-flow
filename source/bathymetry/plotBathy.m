function ph = plotBathy(varargin)
%PLOTBATHY Plots bathymetry into a set of axes
% 
% Syntax:  
%
%       plotBathy() prompts user to select a bathymetry file, loads the
%       file, and draws a patch object containing a triangulation of the seabed
%       elevation / bathymetry coloured according to vertical elevation into the
%       current axes.
%
%       plotBathy(TR) draws the bathymetry contained in the TriRep object
%       TR containing the triangulation of the bathymetry.
%
%       plotBathy(fileName) loads a TriRep object from file pointed to by
%       fileName, which can be an OAS *.tres (with optional bathymetry saved) or
%       a *.mat file containing the object TR.
%
%       plotBathy(AX, ...) plots into the axes AX.
%
%       h = plotBathy(...) returns a handle to the drawn patch object.
%
% Inputs:
%    
%       TR                  TriRep      A TriRep object of any size manageable
%                                       by the current compouter system.
%
%       fileName            string      A filename, including full or relative
%                                       path and extension, of a *.tres or *.mat
%                                       file containing the variable TR, from
%                                       which the TR input is loaded.
%
%       AX                  handle      A handle to a set of axes into which
%                                       bathymetry is plotted (not necessarily
%                                       the current axes).
%
% Outputs:                  
%
%       h                   handle      Handle to the patch object drawn
%
% References:               none
% Future Improvements:      none
% Other m-files required:   none
% Subfunctions:             none
% Nested functions:         none
% MAT-files required:       none
%
% Author:                   T. H. Clark
% Work address:             3 Charles Babbage Road
%                           Cambridge
%                           CB3 0GT
% Email:                    tom.clark@oceanarraysystems.com
% Website:                  www.oceanarraysystems.com
%
% Revision History:         04 August 2014      Created
%                           07 August 2014      Added a missing patch handle
%                                               output

% Parse the inputs
plotHandle = [];
switch nargin
    case 0
        TR = loadBathy();
    case 1
        if ishandle(varargin{1})
            TR = loadBathy();
            plotHandle = varargin{1};
        elseif strcmpi(class(varargin{1}),'TriRep')
            TR = varargin{1};
        else
            TR = loadBathy(varargin{1});
        end
    case 2
        if strcmpi(class(varargin{2}),'TriRep')
            TR = varargin{2};
        else
            TR = loadBathy(varargin{2});
        end
        plotHandle = varargin{1};
end
        
      
if isempty(plotHandle)
    
    ph = trisurf(TR.Triangulation,TR.X(:,1),TR.X(:,2),TR.X(:,3));
    set(ph,'EdgeColor','none')
    axis equal
    
else
    
    % Make the desired axes current
    axes(plotHandle)
    ph = trisurf(TR.Triangulation,TR.X(:,1),TR.X(:,2),TR.X(:,3));
	set(ph,'EdgeColor','none')

end