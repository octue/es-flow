function convertToV73(file)
%CONVERTTOV73 Converts a mat file or files to version 7.3
% Coverting the mat file format allows retrieval of partial variables from a mat
% file as opposed to the requirement of loading all at once in versions pre-7.3.
%
% Syntax:
%       convertToV73(fileName) If fileName is a string, CONVERTTOV73 converts
%       the mat file pointed to by the full or relative fileName to version 7.3.
%       If fileName is a cell array of strings then each of the files needs to
%       be converted. 
%
%       [...] = FunctionName(...,...) Does something else
%
%
% Inputs:
%
%       fileName        string or cell array of strings
%                                   Contains full or relative filename
%
% Outputs:
%
%       output1         [p x q]     Description of output1 including units if 
%                                   any
%
%       outputN         [m x n]     Description of outputN including units if 
%                                   any
%
% See Also: SAVE, LOAD, MATFILE, MATFILES.
%
% References:
%
%   [1] Full references go here with enough detail to allow citation. Ensure the
%       reference document is stored in OAS' Mendeley database so it can be
%       accessed by the whole company.
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
% Revision History:        	21 February 2015  	Created
%                           02 April 2015       Modified to treat all input
%                                               files as .mat regardless of file
%                                               extension. This allows use with
%                                               e.g. OAS' .adcp file format,
%                                               which is actually a mat file.
%                           03 April 2015       Fixed bug in display of the
%                                               number of files to be converted
%                                               (the progress indicator).
%
% Copyright (c) 2014-2015 Ocean Array Systems, All Rights Reserved.


% Handle string/cell inputs
if ischar(file)
    fileCell = {file};
elseif iscell(file)
    fileCell = file;
else
    error('OASUtils:InvalidInputType','Input must be a string or cell array of strings')
end
    
% Convert all the files in the cell
for i = 1:numel(fileCell)
    [~, name] = fileparts(fileCell{i});
    disp(['Converting file: ' name ' (file ' num2str(i) ' of ' num2str(numel(fileCell)) ')'])
    s = load(fileCell{i},'-mat'); %#ok<NASGU>
    save(fileCell{i},'-struct','s','-v7.3')
end


end
