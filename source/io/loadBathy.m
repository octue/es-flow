function [TR] = loadBathy(varargin)
%LOADBATHY Loads bathymetry data into current workspace
% Able to load a triangulation of bathymetry from an OAS tidal resource model
% file (*.tres), an OAS bathymetry file (*.bathy), or a *.mat file containing a
% TriRep object.
%
% Syntax:  
%       [TR] = loadBathy()
%           Select a file containing a triangulation of seabed bathymetry and
%           load the TriRep object into the workspace. Suitable fie types
%           include OAS Tidal RESource files (*.tres), OAS Bathymetry files
%           (*.bathy) and MATLAB (*.mat) data files containing a TriRep object.
%
%       [TR] = loadBathy(fileName)
%           Load the triangulation in fileName into the current workspace
%
% Inputs:
%
%       fileName        string          String containing either the name of a
%                                       blade file on the MATLAB search path, or a
%                                       full (non-relative) path and file name,
%                                       e.g.:
%                                           '/Users/thc29/Documents/myData.tres'
%                                       
% Outputs:
%
%       TR              TriRep          A TriRep object of any size manageable
%                                       by the current computer system,
%                                       containing a triangulation of the
%                                       seabed. The TriRep is not necessarily
%                                       convex or closed.
%
% Future Improvements:      
%
%   [1] Addition of ABPmer standard file formats.
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
% Created:          04 August 2014
% Revisions:        


% Check inputs for validity
if nargin > 2
    error('MATLAB:loadBathy','Too many input arguments')
end

if nargin > 0
    if ischar(varargin{1})
        % Single string
        if (exist(varargin{1},'file') ~= 2) && (exist([varargin{1} '.tres'],'file') ~= 2) && (exist([varargin{1} '.bathy'],'file') ~= 2)  && (exist([varargin{1} '.mat'],'file') ~= 2)
            error('MATLAB:loadBathy',['File not found: ' varargin{1}])
        else
            fileName = varargin{1};
        end
        
    else
        % The input was something other than a string
        error('MATLAB:loadBathy','If used, the first input argument must a character string containing a bathymetry (*.tres, *.bathy or *.mat) data file name')
    end
    
else
    % Use a GUI to determine the file to select
    default_path = retrieve_path('loadBathy');
    
    % Navigate to a different directory, and select an arbitrary file
    [name, path] = uigetfile({'*.tres','OAS TidalRESource files (*.tres)';'*.mat','MAT-files (*.mat)'},'Select a Tidal RESource data file (*.tres)',default_path);

    % Update the database for the most recent path
    update_default_path(path,'loadBathy');
    
    % Set the full filename
    fileName = [path name];
    
end


% Actually load the data
resourceData = load(fileName,'-mat','TR');
TR = resourceData.TR;

% % Display a list of fields
% disp('____________________________')
% disp('loadBathy.m: AVAILABLE FIELDS')
% disp(fields(resourceData));
% disp('____________________________')

    
    
    
    