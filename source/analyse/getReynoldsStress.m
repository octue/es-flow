function [reynoldsStress] = getReynoldsStress(adcpDataWindow, uTau, varargin)
%GETREYNOLDSSTRESS determine Reynolds Stress profile(s) for the current ADCP
% data window. Recommended windows are of order ten minutes.
% Reynolds Stress is:
%     overbar(u_i' u_j')/(uTau^2)
% with i = 1:3, j = 1:3 i.e. there are 9 terms in this stress tensor.
%
% Syntax:  
%       [reynoldsStress uBar] = getReynoldsStress(adcpDataWindow, U_tau)
%           Determines the Reynolds Stress profiles averaged over the time
%           spanned by the data in adcpDataWindow and normalised using the
%           wall friction velocity U_tau.
%
% Inputs:
%
%       adcpDataWindow  structure       An OAS standard adcp data structure,
%                                       windowed to the time period over which
%                                       to average the Reynolds Stresses.
%                                       See help('loadADCP') for fields
%                                       reference.
%
%       uTau            [1 x 1]         Skin friction velocity used to normalise
%                                       the Reynolds stresses 
%
% Outputs:
%
%   	reynoldsStress      [nBins x 6] Reynolds stress profile at nBins
%                                       locations for each of the 9 elements of
%                                       the stress tensor. Note that three terms 
%                                       are redundant in the symmetric tensor
%                                       thus not computed. Columns:
%                                       [r11 r12 r13 r22 r23 r33]
%
%       z                   [nBins x 1] The vertical locations
%
%       timeRange           [2 x 1]     Timestamps (Matlab datenum format) of
%                                       the first and last elements
%
%       nTimes              [1 x 1]     Number of instantaneous velocity
%                                       profiles over which the averaging was
%                                       performed
%
% Future Improvements: 
%
%   [1] Addition of inputs for alternative values of uBar for the Reynolds
%       decomposition
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
% Revisions:        

% Mask out invalid data using any masking flags
if isfield(adcpDataWindow.flags,'invalidMask')
    mask = adcpDataWindow.flags.invalidMask;
    u1 = adcpDataWindow.u(:,~invalidMask);
    u2 = adcpDataWindow.v(:,~invalidMask);
    u3 = adcpDataWindow.w(:,~invalidMask);
else
    u1 = adcpDataWindow.u;
    u2 = adcpDataWindow.v;
    u3 = adcpDataWindow.w;
end

% Perform Reynolds decomposition 
u1Bar = mean(u1,2);
u2Bar = mean(u2,2);
u3Bar = mean(u3,2);
u1Prime = bsxfun(@minus,u1,u1Bar);
u2Prime = bsxfun(@minus,u2,u2Bar);
u3Prime = bsxfun(@minus,u3,u3Bar);

% Calculate the stress tensor
r11 = mean(u1Prime.*u1Prime,2)./(uTau^2);
r12 = mean(u1Prime.*u2Prime,2)./(uTau^2);
r13 = mean(u1Prime.*u3Prime,2)./(uTau^2);
% r21 = mean(u2Prime.*u1Prime,2)./(uTau^2); % redundant
r22 = mean(u2Prime.*u2Prime,2)./(uTau^2);
r23 = mean(u2Prime.*u3Prime,2)./(uTau^2);
% r31 = mean(u3Prime.*u1Prime,2)./(uTau^2); % redundant
% r32 = mean(u3Prime.*u2Prime,2)./(uTau^2); % redundant
r33 = mean(u3Prime.*u3Prime,2)./(uTau^2);

% Concatenate to the output
reynoldsStress = [r11 r12 r13 r22 r23 r33];