function [path] = retrieve_path(varargin)
%RETRIEVE_PATH retrieves saved default pathways.
%   Default directory paths are saved in the database file:
%           matlabroot\default_paths.mat
%   retrieve_path accesses the database and returns either the last used 
%   path, or the last used path of a specific type.
%   If the database file does not exist, a default path (the matlab current
%   directory) is returned. Use update_default_path(), which will create
%   the database file.
%
% Syntax:
%   [path] = retrieve_path
%   [path] = retrieve_path(type)
%
% Inputs:
%
%       type        string      (OPTIONAL)
%                               Contains the type of path to be retrieved. 
%                               Eg if there is a particular directory where
%                               you normally keep your image files, you
%                               could use type = 'imfiles' to save and
%                               retrieve that particular path. If no type
%                               is specified, then the default type,
%                               'last', is used to return the last used
%                               path of no specified type. Thus the only
%                               condition is that you can't use 'last' as
%                               your own type!
%                               N.B. Type strings are case-insensitive, so
%                               'ImFiles' is the same as 'imfiles'
%
% Outputs:
%       path        string      Contains the path retrieved from the 
%                               database. The path is terminated in a
%                               slash, denoting that it is a directory path
%                               rather than a filename.
%                               If the database file could not be found,
%                               or upon first execution of these codes,
%                               then the matlab current directory is
%                               returned.
%
% Examples:
%       % See update_default_path.m for example code
%       help('update_default_path')
%       
% Other m-files required:   update_default_path.m
% Subfunctions:             none
% Nested Functions:         none
% MAT-files required:       /default_paths.mat 
%                               (automatically created)
%
% Author:           T.H. Clark
% Work address:     Fluids Laboratory
%                   Cambridge University Engineering Department
%                   2 Trumpington Street
%                   Cambridge
%                   CB21PZ
% Email:            t.clark@cantab.net
% Website:          http://cambridge.academia.edu/ThomasClark/
%
% Created:          07 August 2009 
% Last revised:     19 August 2009

% Error handling to prevent disruption of your code
rootdir = which('update_default_path');
rootdir = rootdir(1:end-21);

try
    % Load the default paths from file
    def_dirs = load([rootdir 'default_paths.mat']);
    
    % Get the required type
    if nargin > 0;      type = lower(varargin{1});
    else                type = 'last';
    end
    
    % Get the path using dynamic fieldname...
    path = def_dirs.(type);
    
    % If the last character of the path string is not \ or / then add one
    % so that functions like uigetfile() recognise the path as being to a
    % directory, rather than a file
    if ~strcmp(path(end),'\') && ~strcmp(path(end),'/')
        path = [path filesep];
    end
catch %#ok<CTCH>
    % Prevent any fuss if the load is unsuccessful. It probably means that
    % update_default_directory has not yet been run, therefore the database
    % file isn't initialised yet.
    path = cd;
end




