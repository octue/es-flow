function [binInds] = structHist(data, varargin)
%STRUCTHIST Sorts data structures into bins based on one or more scalar properties
%
% Syntax:
%       [bins, binInds] = binWindows(data, 'Property1', binEdges1)
%
%       [bins, binInds] = binWindows(..., 'PropertyN', binEdgesN)
%
%
% Inputs:
%
%       data         [any]       A multidimensional structure with nWin
%                                   entries. Each element of the data structure 
%
%       'PropertyN'     string      Fieldname of a scalar property of data,
%                                   used to sort different elements of in
%                                   data into different bins.
%
%       binEdgesN       [1 x nBins+1] 
%                                   Edges of bins for each property.
%
% Outputs:
%
%       bins            [nBins1 x ... x nBinsN]
%                                   Multidimensional structure containing an
%                                   element for each bin 
%
%       binInds         [nWin x nProps]
%                                    The indices relating each window to the bin
%                                    in which it resides for each property
%
% See Also: HIST
%
% Future Improvements:
%
%   [1] Disgusting code. Tidy up!
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
% Revision History:        	25 February 2015    Created
%                           27 February 2015    Documented
%                           11 March 2015       Fixed a bug where index could
%                                               exceed the number of bins in the
%                                               input edges vector.
%
% Copyright (c) 2014-2015 Ocean Array Systems, All Rights Reserved.

% Sort out the bins
binVariables = varargin(1:2:end);
binEdges = varargin(2:2:end);

% Check that bin edges are provided as strictly increasing with no repetition
for iVar = 1:numel(binEdges)
    for iBin = 1:numel(binEdges)
        if ~issorted(binEdges{iBin}) || (numel(unique(binEdges{iBin})) ~= numel(binEdges{iBin}))
            error('binEdges must be strictly increasing with no repeating values')
        end
    end
end
    
% Preallocate variables to hold bin values
binValues = zeros(size(binVariables));

% Assign bin indices array for each window
binInds = zeros(numel(data),numel(binVariables));

% Initialise progress indicator
dispstat('','init')
dispstat('structHist.m: Beginning binning process','timestamp')
for iWin = 1:numel(data)
    dispstat(['Binning structure entry ' num2str(iWin) ' of ' num2str(numel(data))],'timestamp')
    for iVar = 1:numel(binVariables)
        binValue = data(iWin).(binVariables{iVar});
        if isempty(binValue)
            warning(['Window ' num2str(iWin) ': Empty entries found in the sorting fields of input data structure. Assigning NaN value to bin index.'])
            binInds(iWin,iVar) = NaN;
        elseif numel(binValue) ~= 1
            error('Attempt to bin on vector values. You must bin based on single value fields in the input data structure.')
        elseif imag(binValue) ~= 0
            warning('Complex values in sort fields detected. Assigning NaN value to bin index.')
            binInds(iWin,iVar) = NaN;
        else
            binValues(iVar) = binValue;
            ind = find(binValue >= binEdges{iVar}, 1, 'last');
            if isempty(ind) || (ind >= numel(binEdges{iVar}))
                warning(['Bin range excluding window ' num2str(iWin) ' with ' binVariables{iVar} ' = ' num2str(binValue)])
            else
                binInds(iWin,iVar) = ind;
            end
        end
    end

end


end % End main function