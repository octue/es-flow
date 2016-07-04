function [normals, slope] = bathyGradients(TR)
%BATHYGRADIENTS Returns normal directions and slope (in degrees from vertical)
%for a 2D triangulation. Useful for checking the slope of a seabed in
%preparation for planning an ADCP deployment.
% 
% Syntax:  
%
%       [normals] = bathyGradients(TR) Returns the unit normal vector to
%       triangles within the 2D triangulation TR. Trivially overloads the
%       faceNorm method of MATLAB's triangulation class.
%
%       [normals, slope] = triBathy(...) Returns maximum slope (in degrees from
%       the upward vertical) of each face in TR.
%
% Inputs:
%
%       TR                  [1 x 1] TriRep Object
%                                               Triangulated surface of the
%                                               bathymetry data. See
%                                               help('TriRep') for more details.
%
% Outputs:
%
%       normals             [nT x 3]            Components (X,Y,Z) of the unit
%                                               normal vector for each 2D
%                                               simplex within the TR
%                                               triangulation 
%
%       slope               [nT x1]             Angle in degrees between the
%                                               unit normal vector and the
%                                               vertical [0 0 1], for each
%                                               simplex in triangulation TR.
%
% References:               none
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
% Revision History:         06 August 2014      Created


% Compute surface normals of the triangulation in Tri, trivial overload of the
% faceNorm method:
normals = faceNormals(TR);

% Compute the angle between this vector and the upward vertical
d = sqrt(sum(normals(:,1:2).^2,2));
slope = atand(d./normals(:,3));

