function [adcpData] = loadADCP(varargin)
%loadADCP Loads ADCP data into the current workspace
% Able to load seabed mounted ADCP data from OAS format ADCP data files (*.adcp)
%
% Syntax:  
%       [adcpData] = loadADCP()
%           Select a file containing ADCP data and load it into the MATLAB
%           workspace.
%
%       [adcpData] = loadADCP(fileName)
%           Load the ADCP data in fileName into the current workspace
%
% Inputs:
%
%       fileName        string          String containing either the name of a
%                                       blade file on the MATLAB search path, or a
%                                       full (non-relative) path and file name,
%                                       e.g.:
%                                           '/Users/thc29/Documents/myData.mat'
%                                       
% Outputs:
%
%   adcpData         	 structure      A structure containing velocity data in
%                                       the following form:
%
%             .u         [nBins x nT]   Velocity for nBins bins above
%                                       the seabed mounted ADCP for nT time
%                                       points. Units m/s. Default direction is
%                                       Eastings velocity.
%
%             .v         [nBins x nT]   Velocity as per u. Units m/s. Default
%                                       direction is Northings velocity
%
%             .w         [nBins x nT]   Vertical velocity as per u, positive
%                                       upward, making a right handed set with u
%                                       and v. Units m/s.
%
%             .z         [nBins x 1]    Bin heights above the seabed Units m
%
%             .zUnit     [1 x 1]        ADCP device height above seabed in m
%
%             .beamAngle [1 x 1]        Beam angle in degrees, typically 20 or
%                                       25.
%             
%             .direction [1 x 1]        Orientation of the u direction, in
%                                       degrees anticlockwise from East. Default
%                           or          0 gives u as Eastings velocity, v as
%                                       Northings velocity. Direction can change
%                        [1 x nT]       over time (specified for nT times),
%                                       typically allowing the primary direction
%                                       to be aligned always with the bulk flow
%                                       (as determined by a moving average or
%                                       similar filter - see flowDirection.m).
%
%             .d         [1 x nT]       Water depth (from the seabed) at nT
%                                       times. Units m.
%
%             .t         [1 x nT]       Time series containing the datenum of
%                                       acquisition. Times must be strictly
%                                       increasing but need not be monotonic.
%       
%             .flags     structure      Contains notes/warning flags and text
%                                       descriptions of their meanings
%
% Future Improvements:      
%
%   [1] Addition of ADCP manufacturer ASCII file formats (e.g. as written by RDI
%       Workhorse)
%
%   [2] Auto recognition of file types 
%
%   [3] Data types using multiple files (through a fileName cell, or through
%       autorecognition of the header file type
%
%   [4] Limitation to a timestamp range and/or an index range
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
% Created:          07 July 2014
% Revisions:        23 July 2014        Added beamAngle and zUnit allowing the
%                                       spectrumADCP function to compute the
%                                       beam separation frequency and scales
%                   13 August 2014      Changed text output to remove list of
%                                       fields but display a 'Loaded <filename>'
%                                       message


% Check inputs for validity
if nargin > 2
    error('MATLAB:loadADCP','Too many input arguments')
end

if nargin > 0
    if ischar(varargin{1})
        % Single string
        if (exist(varargin{1},'file') ~= 2) && (exist([varargin{1} '.adcp'],'file') ~= 2) && (exist([varargin{1} '.mat'],'file') ~= 2)
            error('MATLAB:loadADCP',['File not found: ' varargin{1}])
        else
            fileName = varargin{1};
        end
        
    else
        % The input was something other than a string
        error('MATLAB:loadADCP','If used, the first input argument must a character string containing an ADCP (*.adcp or *.mat) data file name')
    end
    
    if nargin > 1
        if ischar(varargin{2})
            % String specifying file type
            fileType = varargin{2};
        else
            % The input was something other than a string
            error('MATLAB:loadADCP','If used, the second input argument must a character string denoting the type of file being read. Valid strings are ''partrac'' ''oas'' (more to follow). Case insensitive.')
        end
    else
        % Default to OAS file type
        fileType = 'oas';
    end
    
else
    % Use a GUI to determine the file to select
    default_path = retrieve_path('loadADCP');
    
    % Navigate to a different directory, and select an arbitrary file
    [name, path] = uigetfile('*.adcp','Select an ADCP data file (*.adcp)',default_path);

    % Update the database for the most recent path
    update_default_path(path,'loadADCP');
    
    % Set the full filename
    fileName = [path name];
    
end


% Actually load the data
switch lower(fileType)
    case 'oas'
        if ~strcmp(fileName(end-4:end),'.adcp') 
            fileName = [fileName '.adcp'];
        end
    case 'partrac'
        if ~strcmp(fileName(end-3:end),'.mat') 
            fileName = [fileName '.mat'];
        end
    otherwise
        error('MATLAB:InvalidInput','File Type input does not describe a recognised file type. Try ''oas'' or ''partrac''')
end
adcpData = load(fileName,'-mat');

% % Display a list of fields
% disp('____________________________')
% disp('loadACDP.m: AVAILABLE FIELDS')
% disp(fields(adcpData));
% disp('____________________________')

% We don't strip the structure of fields that are not mandated, to avoid
% inadvertent data loss. However, we do error if required fields are not present
requiredFields = {'u';'v';'w';'z';'zUnit';'beamAngle';'direction';'d';'t';'flags'};
for iCtr = 1:numel(requiredFields)
    if ~any(isfield(adcpData,requiredFields{iCtr}))
        error('MATLAB:InvalidFile','ADCP data file is missing one or more of the following fields: ''u'' ''v'' ''w'' ''z'' ''zUnit'' ''beamAngle'' ''direction'' ''d'' ''t'' ''flags''')
    end
end

% Display progress
disp(['loadADCP: Loaded file ' fileName])

