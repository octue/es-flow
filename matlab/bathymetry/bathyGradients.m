function [normals, slope] = bathyGradients(TR, flag)
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
%       [normals, slope] = bathyGradients(TR) Returns maximum slope (in degrees
%       from the upward vertical) of each face in TR.

%       [...] = bathyGradients(TR, flag) Returns maximum slope (in degrees
%       from the upward vertical) of each face in TR if flag = 'faces', or for
%       each vertex in TR if flag = 'vertices'
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
% Author:                   T. H. Clark
% Work address:             Ocean Array Systems Ltd
%                           Hauser Forum
%                           3 Charles Babbage Road
%                           Cambridge
%                           CB3 0GT
% Email:                    tom.clark@oceanarraysystems.com
% Website:                  www.oceanarraysystems.com
%
% Copyright (c) 2016 Ocean Array Systems, All Rights Reserved.

% COnvert the triangulation type if it's the outdated TriRep
if isa(TR, 'TriRep')
    TR = triangulation(TR.Triangulation,TR.X);
end

% Set default behaviour to deliver faces
if nargin == 1
    flag = 'faces';
end

% Get normal vectors either at faces or at vertices
switch lower(flag)
    
    case 'faces'
        % Compute surface normals of the triangulation in Tri
        normals = faceNormal(TR);

    case 'vertices'
        % Compute surface normals of the triangulation in Tri at vertices
        normals = vertexNormal(TR);
        
    otherwise
        error('Unrecognised flag string. Try ''face'' or ''vertices''.')
        
end

% Compute the angle between these vectors and the upward vertical
d = sqrt(sum(normals(:,1:2).^2,2));
slope = atand(d./normals(:,3));

end


