function [resourceData] = loadResource(varargin)
%loadResource Loads Tidal Resource model into current workspace
% Able to load 2D resource model data (e.g. produced by Mike21) from OAS format
% tidal resource data files (*.tres).
%
% Syntax:  
%       [resourceData] = loadResource()
%           Select an OAS Tidal RESource (*.tres) file using a GUI and load
%           model data into a structure.
%
%       [resourceData] = loadResource(fileName)
%           Load the resource data in fileName.
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
%   resourceData     	structure       A structure containing velocity data in
%                                       the following form:
%
%       .x              [nY x nX]       Regular spatial grid (in meshgrid form)
%                                       on which velocity values are calculated.
%                                       X contains longitude values.
%
%       .y              [nY x nX]       Latitude as per field X.
%
%       .z              [nY x nY]       Height of the sea bed in m relative to
%                                       some datum (see 'zDatum'). Positive
%                                       upward.
%
%       .zDatum         string          Indicates the vertical datum used for
%                                       bathymetry (Z). Typically 'MSL' (Mean
%                                       Sea Level).
%
%       .gridUnits      string          String containing the coordinate system
%                                       in which X and Y are given. 'utm29N' is
%                                       typical for UK regions.
%
%   	.u              [nY x nX x nT]  Eastings velocity in m/s on a regular
%                                       grid described by fields X, Y, for nT
%                                       points in time described by field t.
%
%   	.v              [nY x nX x nT]  Northings velocity as per u.
%
%       .depth          [nY x nX x nT]  Water depths in m as per u.
%
%       .t              [1 x nT]        MATLAB standard for datenum describing
%                                       the time for which each step of the
%                                       model applies.
%
%     	.tInds          [1 x nT]        (Optional)
%                                       Where a subset (in time) of a wider
%                                       model is stored, these indices refer to
%                                       the entries in the wider model from
%                                       which the present dataset was extracted.
%
%       .bathymetry     TriRep          (Optional)
%                                       A MATLAB TriRep object containing the
%                                       site bathymetry for the gridded region.
%
% Future Improvements:      none
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
% Created:          01 August 2014
% Revisions:        04 August 2014      Added depth, Z and zDatum fields and
%                                       tidied redundant filetype capability in
%                                       the code.


% Check inputs for validity
if nargin > 1
    error('MATLAB:loadResource','Too many input arguments')
end

if nargin > 0
    if ischar(varargin{1})
        % Single string
        if (exist(varargin{1},'file') ~= 2) && (exist([varargin{1} '.tres'],'file') ~= 2) && (exist([varargin{1} '.mat'],'file') ~= 2)
            error('MATLAB:loadResource',['File not found: ' varargin{1}])
        else
            fileName = varargin{1};
        end
        
    else
        % The input was something other than a string
        error('MATLAB:loadResource','If used, the first input argument must a character string containing a resource (*.tres or *.mat) data file name')
    end
    
else
    % Use a GUI to determine the file to select
    default_path = retrieve_path('loadResource');
    
    % Navigate to a different directory, and select an arbitrary file
    [name, path] = uigetfile({'*.tres','OAS TidalRESource files (*.tres)';'*.mat','MAT-files (*.mat)'},'Select a Tidal Resource model file', default_path);

    % Update the database for the most recent path
    update_default_path(path,'loadResource');
    
    % Set the full filename
    fileName = fullfile(path,name);
    
end

% Actually load the data
resourceData = load(fileName,'-mat');
