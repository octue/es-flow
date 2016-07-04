function [adcpWindow] = conditionalAverage(adcpData, bin, quadrant)
%CONDITIONAVERAGE
%
%
% Syntax:  
%       [adcpWindow] = getWindow(adcpInput, mode, window)
%           
%
% Inputs:
%
%       adcpInput       structure       An OAS standard adcp data structure.
%                                       See help('loadADCP') for fields
%                                       reference.
%
%       quadrant        
%
%       bin          
%
%
% Outputs:
%
%   	adcpWindow      structure       An OAS standard adcp data structure
%                                       extracted from a limited part of the
%                                       input as specified by the window limits
%                                       or indices.
%                                       See help('loadADCP') for fields
%                                       reference.
%
%
% Future Improvements: 
%
%   [1] Handling of additional or user customised fields
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
% Created:          10 July 2014
% Revisions:        


% Get the indices to extract. NB don't use a logical mask, since we wish to
% allow reordering through the indices vector.
switch lower(mode)
    case 'indexrange'
        
        % Handle noninteger inputs and inputs of an incorrect size
        if any(fix(window) ~= window) || (numel(window) ~= 2)
            error('MATLAB:InvalidInput','if using the ''indexRange'' mode, the window must be specified as a pair of integer indices')
        end
        
        % Get linear index series
        indices = min(window):1:max(window);
        
    case 'indices'
        
        % Handle noninteger inputs
        if any(fix(window) ~= window)
            error('MATLAB:InvalidInput','if using the ''indices'' mode, the window must be specified as series of integer indices')
        end
        
        % Use the input directly as index series
        indices = window;
        
    case 'timerange'
        
        % Get indices corresponding to entries within the desired time range
        indices = find((adcpInput.t >= min(window)) && (adcpInput.t <= max(window)));
        
    otherwise
        error('MATLAB:InvalidInput','Invalid mode string input. Try ''indexRange'', ''indices'' or ''timeRange'' (not case sensitive)')
        
end


% Handle the eventuality of an empty window
if isempty(indices)
    error('MATLAB:InvalidInput','Window contains no entries. Try widening the index or time range selected')
end


% Create new structure with windowed fields
adcpWindow.u = adcpInput.u(:,indices);
adcpWindow.v = adcpInput.v(:,indices);
adcpWindow.w = adcpInput.w(:,indices);
adcpWindow.d = adcpInput.d(:,indices);
adcpWindow.t = adcpInput.t(:,indices);

% Add the fields that don't get windowed
adcpWindow.z = adcpInput.z;

% We need to find flags that have dimension [anything x nT] and window them too,
% otherwise just copy the whole flag field
flagFields = fields(flags);
for i = 1:numel(flagFields)
    if size(adcpInput.flags.(flagFields{i}),2 == nT)
        adcpWindow.flags.(flagFields{i}) = adcpInput.flags.(flagFields{i})(:,indices);
    else
        adcpWindow.flags.(flagFields{i}) = adcpInput.flags.(flagFields{i});
    end
end
        

end % End main function