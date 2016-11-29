function [TR] = triBathy(dataFile, varargin)
%TRIBATHY Reads CSV data and triangulates bathymetry data. 
% Outputs the triangulation in a MATAB TriRep object; optionally saves outputs
% to a .STL file and/or a .MAT file.
% 
% Syntax:  
%
%       TR = triBathy(dataFile)
%   
%            Uses point depth locations contained in a csv file (e.g. as
%            downloadable from Ref [1]) to create a triangulation of seabed
%            surface
%
%       TR = triBathy(dataFile, 'Parameter', value, ...)
%
%            Uses non-default options for processing, scaling and saving the
%            input data.
%
% Inputs:
%    
%       dataFile             string     String containing full or relative path
%                                       (including file extension) of the .CSV
%                                       or similar file to read in. Number of
%                                       header lines is automatically detected.
%                                       File must contain 3 columns (Eastings,
%                                       Northings and Height translating to
%                                       X,Y,Z locations of points on the seabed)
%
% Optional Inputs:
%
%       Parameter           Value
%       'Delimiter'         string      Delimiter to use in reading the file,
%                                       strings as per importData(). Leaving the
%                                       string empty (default) attempts to
%                                       determine the delimiter automatically.
%
%       'NHeaderLines'      [1 x 1]     Only used if 'Delimiter' also specified.
%                                       Number of header lines to be skipped.
%                                       Leaving 'NHeaderLines' or the
%                                       'Delimiter' parameter unspecified
%                                       results in attempt to auto-detect the
%                                       number of header lines.
%
%       Scale               [1 x 3] double      Default [1 1 1]
%                                               Scale factor applied to each
%                                               dimension of the input data
%
%       OriginOffset        [1 x 3] double      Default [0 0 0]
%                                               Offsets the origins of the data
%                                               in 3 dimensions
%
%       OriginOffsetUnits   string              Default 'input'
%                                               Dictates the units in which the
%                                               origin offsets are given. Valid
%                                               strings are 'input' and 'output'
%
%       Alpha               [1 x 1] double      Default [Inf]
%                                               Specifying a value for alpha
%                                               invokes an alpha shape algorithm
%                                               to trim edges of the dataset.
%                                               This is useful for non-convex
%                                               regions of bathymetry data, in
%                                               which non-convexities are filled
%                                               in by the default triangulation
%                                               algorithm (delaunay). 
%                                               Alpha specifies the radius of a
%                                               hypothetical circle which, when
%                                               rolled around the outside of the
%                                               region, would be prevented from
%                                               entering. Good typical values
%                                               are of order 10x - 20x the
%                                               sampling distance of the data.
%
%       AlphaUnits          string              Default 'input'
%                                               Specifies whether the Alpha
%                                               radius value is given in the
%                                               same units as the input data or
%                                               the scaled output data.
%
%       SaveSTL             [1 x 1] logical     Default false
%                           string              If a character string is given,
%                                               STL file is saved to the
%                                               location and filename given in
%                                               the string.
%
%       SaveMAT             [1 x 1] logical     Default false
%                              or               Same behaviour as SaveSTL, but
%                           string              saves a TriRep object to a
%                                               MATLAB .MAT binary file type.
%   
%       XLim                [1 x 2] double      [-Inf Inf];
%       YLim                [1 x 2] double      [-Inf Inf];
%                                               Apply limits to the dataset to
%                                               triangulate; i.e. discard data
%                                               points outside the maxima and
%                                               minima in these X and Y
%                                               directions.
%
%       LimUnits            string              Default 'input'
%                                               Specify units for XLim and YLim
%                                               parameters in the  same way as
%                                               AlphaUnits above.
%
% Outputs:
%
%       TR                  [1 x 1] TriRep Object
%                                               Triangulated surface of the
%                                               bathymetry data. See
%                                               help('TriRep') for more details.
%
% References:
%
%   [1] UK Hydrographic Office
%
% Future Improvements:      none
% Other m-files required:   none
% Subfunctions:             none
% Nested functions:         none
% MAT-files required:       none
%
%
% Author:                   T. H. Clark
% Work address:             3 Charles Babbage Road
%                           Cambridge
%                           CB3 0GT
% Email:                    tom.clark@oceanarraysystems.com
% Website:                  www.oceanarraysystems.com
%
% Revision History:         17 March 2014       Created
%                           18 March 2014       Input options refined and
%                                               documented
%                           04 August 2014      Removed plot options and
%                                               associated code which is now
%                                               easily accessed through
%                                               replacement function plotBathy()

% Set and behaviour options using input pairs
opts.Delimiter = [];
opts.NHeaderLines = [];
opts.Scale = [1 1 1];
opts.OriginOffset = [0 0 0];
opts.OriginOffsetUnits = 'input'; % other option is output
opts.Alpha = Inf;
opts.AlphaUnits = 'output';
opts.SaveSTL = false;
opts.SaveMAT = false;
opts.XLim = [-Inf Inf];
opts.YLim = [-Inf Inf];
opts.LimUnits = 'output';
if nargin>1
    opts = parse_pv_pairs(opts,varargin);
end

% Read the CSV to a variable. Use importData as it automatically handles header
% lines.
if isempty(opts.Delimiter)
    [X, delimiterOut, headerLinesOut] = importdata(dataFile);
    disp('triBathy.m: Automatically detected delimiter and header lines')
elseif isempty(opts.NHeaderLines)
    [X, delimiterOut, headerLinesOut] = importdata(dataFile, opts.Delimiter);
else
    [X, delimiterOut, headerLinesOut] = importdata(dataFile, opts.Delimiter, opts.NHeaderLines);
end

% Multiple output cases from importdata... handle:
if isstruct(X)
    X = X.data;
end 

% Display result:
disp(['triBathy.m: Successfully read file ' dataFile])
disp(['triBathy.m: Data read is of size ' num2str(size(X))])
disp(['triBathy.m: Detected delimiter of type ' delimiterOut])
disp(['triBathy.m: Skipped ' num2str(headerLinesOut) ' header lines.'])

% Display error if data is not 3 columns (x,y,depth)
if size(X,2) ~= 3
    error('MATLAB:triBathy','triBathy.m: Expected data with 3 columns (x, y, elevation)')
end


% Perform unit conversions to output
if strcmpi(opts.OriginOffsetUnits,'input')
    opts.OriginOffset = opts.OriginOffset*opts.Scale(1);
end
if strcmpi(opts.AlphaUnits,'input')
    opts.Alpha = opts.Alpha*opts.Scale(1);
end
if strcmpi(opts.LimUnits,'input')
    opts.XLim = opts.XLim*opts.Scale(1);
    opts.YLim = opts.YLim*opts.Scale(2);
end

% Scale and offset
X = bsxfun(@times, X, opts.Scale);
X = bsxfun(@plus, X, opts.OriginOffset);

% Remove parts of dataset outside limits
Xbool = X(:,1)<min(opts.XLim) | X(:,1)>max(opts.XLim);
Ybool = X(:,2)<min(opts.YLim) | X(:,2)>max(opts.YLim);
X(Xbool|Ybool,:) = [];

% Alpha around the edges if requested
if opts.Alpha ~= Inf;
    
    % Extract the alpha shape
    [~, S] = alphavol(X(:,1:2),opts.Alpha);
    
    % Triangulate the surface with constrained boundary of the alpha shape
    DT.Triangulation = S.tri;
    
else
    
    % Triangulate the surface with unconstrained boundary
    DT = DelaunayTri(X(:,1:2));
    
end

% Output the STL File (binary by default). Allows use of character string to
% define output filename otherwise reverts to default same name as CSV file
if ischar(opts.SaveSTL)
    stlwrite(opts.SaveSTL, DT.Triangulation, X);
end

% Create a TriRep with the 3D, rather than 2D, facets for output
TR = TriRep(DT.Triangulation, X);

% Optionaly save File (v7.3 enabled binary) for further processing
if ischar(opts.SaveMAT)
    save(opts.SaveMAT, 'TR');
end






end

