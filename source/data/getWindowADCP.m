function [adcpWindow, indexRange] = getWindowADCP(adcpInput, mode, window)
%GETWINDOWADCP Gets a window from an adcp data structure
%
% Syntax:  
%       [adcpWindow] = getWindowADCP(adcpInput, mode, window)
%       
% Inputs:
%
%       adcpInput       structure       An OAS standard adcp data structure.
%                                       See help('loadADCP') for fields
%                                       reference.
%
%       mode            string          The mode by which the window is
%                                       described. Possible modes:
%
%                       'indexRange'    Extracts window between two indices of
%                                       entries in the input data structure
%                       'indices'       Extracts window comprising profiles at
%                                       the indices in the input vector
%                       'timeRange'     Extracts window comprising any profiles
%                                       in the adcp data structure taken between
%                                       two times (including profiles at those
%                                       times)
%
%       window          Dependent on mode input:
%                       
%                       [1 x 2]         'indexRange'
%                                       Start and end indices of the data window
%                                       into the input data structure
%
%                       [1 x nP]        'indices'
%                                       Vector of nP indices into the input data
%                                       structure profiles that will comprise
%                                       the window
%       
%                       [1 x 2]         'timeRange'
%                                       Start and end times of the desired
%                                       window in MATLAB datenum format (see
%                                       doc('datenum'). Where profiles in the
%                                       input structure correspond exactly with
%                                       the start and end times, these profiles
%                                       are included.
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
%       indexRange      [1 x 2]         Indices of the first and last entries
%                                       in the input dataset which are included
%                                       within the window. In the case of
%                                       operating in 'indexRange' mode, this is
%                                       specified as an input. Slight
%                                       modification to the way the flag fields
%                                       are indexed to allow ADCP objects to be
%                                       matfile linked.
%
% Future Improvements:      none
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
% Revisions:        11 August 2014      Addition of indexRange output
%                   06 April 2015       Altered documentation and function name
%                                       to reflect the filename change to
%                                       getWindowADCP

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
        indices = find((adcpInput.t >= min(window)) & (adcpInput.t <= max(window)));
        
    otherwise
        error('MATLAB:InvalidInput','Invalid mode string input. Try ''indexRange'', ''indices'' or ''timeRange'' (not case sensitive)')
        
end


% Handle the eventuality of an empty window
if isempty(indices)
    error('MATLAB:InvalidInput','Window contains no entries. Try widening the index or time range selected')
end

% Output the index range
indexRange = [min(indices) max(indices)];

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
flagFields = fields(adcpInput.flags);
nT = numel(adcpInput.t);
flags = adcpInput.flags;
for i = 1:numel(flagFields)
    if size(flags.(flagFields{i}),2) == nT
        adcpWindow.flags.(flagFields{i}) = flags.(flagFields{i})(:,indices);
    else
        adcpWindow.flags.(flagFields{i}) = flags.(flagFields{i});
    end
end
        
% Add the data type field
adcpWindow.DataType = 'ADCP';

end % End main function