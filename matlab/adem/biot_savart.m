%BIOT_SAVART applies a biot savart relation
%
% Inputs:
%     start_nodes                      [3 x N] double
%     end_nodes                        [3 x N] double
%     control_point_locations          [3 x P] double
%     gamma                            [1 x N] double
%     effective_core_radius_squared    [1 x N] double
%     mode                             string, optional, either 'naive' or 'tree'
%
% Outputs:
%     induction [3 x P] double
%
% Author:                   T. Clark
% Work address:             Octue Ltd
%                           Hauser Forum
%                           3 Charles Babbage Road
%                           Cambridge
%                           CB3 0GT
% Email:                    tom@octue.com
% Website:                  www.octue.com
%
% Copyright (c) 2013-8 Octue Ltd. All Rights Reserved.
function induction = biot_savart(startNodes, endNodes, controlPointLocations, gamma, effectiveCoreRadiusSquared, mode)
