function [adcpInput] = rotateADCP(adcpInput, u1Direction)
%ROTATEADCP Rotates ADCP velocities u,v about the vertical to a new coordinate
%frame, typically aligned with the direction of flow.
%
% Syntax:
%
%       [adcpOutput] = rotateADCP(adcpInput, u1Direction)
%       
% Inputs:
%
%       adcpInput       structure       An OAS standard adcp data structure.
%                                       See help('loadADCP') for fields
%                                       reference.
%
%       u1Direction     [1 x 1]         Direction, in degrees anticlockwise from
%                                       East, of the output u1. Can be given as 
%                        or             a single value or as a time series. If a
%                                       time series, nT = numel(adcpInput.t)
%                       [1 x nT]        (i.e. each entry corresponds with a
%                                       profile in adcpInput).
%
% Outputs:
%
%   	adcpOutput      structure       An OAS standard adcp data structure
%                                       of the same size as the input,
%                                       containing additional fields .u1 and .u2
%                                       where 'u1' is aligned with teh input
%                                       direction, and [u1 u2 w] form a
%                                       cartesian, right handed set.
%                                       See help('loadADCP') for fields
%                                       reference.
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
% Revision History:        	15 July 2014        Created
%                           06 April 2015       Modified header to OAS standard
%                                               and changed output so that
%                                               velocities are not duplicated.
%
% Copyright (c) 2014-2015 Ocean Array Systems, All Rights Reserved.

% Standard rotation of a vector within a coordinate frame
% [cos -sin
%  sin  cos]

% To rotate the coordinate frame onto the vector, we do it in the opposite
% direction:
sinDir = sind(-u1Direction);
cosDir = cosd(-u1Direction);

u = cosDir.*adcpInput.u - sinDir.*adcpInput.v;

v = sinDir.*adcpInput.u + cosDir.*adcpInput.v ;

adcpInput.u = u;
adcpInput.v = v;

end % end main function