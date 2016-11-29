function [adcpOutput, mask] = removeNaNs(adcpInput, varargin)
%REMOVENANS Removes rows of ADCP data which contain NaNs
%
% Syntax:
%       [adcpData] = removeNaNs(adcpData) Remove all rows in the input adcpData
%       structure that contain NaN values
%
%       [adcpData] = removeNaNs(adcpData, pc) Removes rows in the input adcpData
%       structure that contain more than pc% NaN values; replaces remaining NaNs
%       with 5 point local median values
%
% Inputs:
%
%       adcpData        structure   ADCP data structure     
%
% Outputs:
%
%       adcpData        structure   ADCP data structure with NaN columns
%                                   removed.
%
%         mask          [nz x 1]    Logical mask, with entries corresponding to 
%                                   bins of the input data. Set to true where
%                                   rows have been removed and false where
%                                   they've been retained.
%
% See Also: correctTime, filterADCP
%
% Future Improvements:
%
%   [1] Generalise to a row removal tool.
%
% References:               none
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
% Revision History:        	06 April 2015       Created
%
% Copyright (c) 2014-2015 Ocean Array Systems, All Rights Reserved.


% Get the indices of columns to delete
if nargin > 1
    p = varargin{1}/100;
else
    p = 0;
end
mask =        sum(isnan(adcpInput.u),2)/size(adcpInput.u,2) > p;
mask = mask | sum(isnan(adcpInput.v),2)/size(adcpInput.v,2) > p;
mask = mask | sum(isnan(adcpInput.w),2)/size(adcpInput.w,2) > p;

% Handle the eventuality of an empty window
if all(mask == true)
    error('MATLAB:InvalidInput','ADCP data contains a NaN in every row')
end

% See future improvements: we'll eventually use indices but for now just use a
% logical mask to give the indices to keep
indices = 1:size(adcpInput.u,1);
indices = indices(~mask);

% Create new structure with windowed fields
adcpOutput.u = adcpInput.u(indices,:);
adcpOutput.v = adcpInput.v(indices,:);
adcpOutput.w = adcpInput.w(indices,:);
adcpOutput.z = adcpInput.z(indices,1);

% Add variables that don't get masked
if isfield(adcpInput,'d')
    adcpOutput.d = adcpInput.d;
end
adcpOutput.t = adcpInput.t;

% We need to find flags that have dimension [n x anything] and window them too,
% otherwise just copy the whole flag field
if isfield(adcpInput,'flags')
    flagFields = fields(adcpInput.flags);
    n = size(adcpInput.u,1);
    flags = adcpInput.flags;
    for i = 1:numel(flagFields)
        if size(flags.(flagFields{i}),1) == n
            adcpOutput.flags.(flagFields{i}) = flags.(flagFields{i})(indices,:);
        else
            adcpOutput.flags.(flagFields{i}) = flags.(flagFields{i});
        end
    end
end

% Retain datatype
if isfield(adcpInput,'DataType')
    adcpOutput.DataType = adcpInput.DataType;
end

% Replace any remaining NaNs with 5 point median values
clearvars adcpInput %saves memory
adcpOutput.u = inpaint_nans(adcpOutput.u);
adcpOutput.v = inpaint_nans(adcpOutput.v);
adcpOutput.w = inpaint_nans(adcpOutput.w);


end
