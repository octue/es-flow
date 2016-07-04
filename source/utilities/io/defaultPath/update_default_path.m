function update_default_path(path,varargin)
%UPDATE_DEFAULT_PATH Creates or updates a file containing default paths
% Creates or updates(in the matlab root directory) a file named
% default_paths.mat, to which is saved a structure containing a number of
% file paths. 
%
% These paths can then be retrieved (see retrieve_path.m) for use with
% functions like uigetfile().
%
% Browse the help for uigetfile() - the option DefaultPath can be used to
% open the file selection dialogue box in a directory other than MATLAB's
% current directory. Each time you use uigetfile or a similar function,
% use update_default_path to save the path of the selected file.
%
% You can associate the path you are saving with a type, too. So, say you
% have two different m-files. One of them uses calibration files (which you
% store in one place) and the other loads image files (which you store in a
% different place). Specifying types allows you to retrieve the most recent
% path used by a specific function
%
% Syntax:  
%       update_default_path(path)
%       update_default_path(path, type)
%
% Inputs:
%       path        string      Contains the path to be saved to the
%                               database. Recommend the use of full,
%                               non-relative pathnames. It doesn't matter
%                               whether the path is terminated in a slash
%                               or not.
%
%       type        string      (OPTIONAL)
%                               Contains the type of path. Eg if you're
%                               saving the path to the directory where you
%                               normally keep your image files, you could
%                               use type = 'imfiles' to save and retrieve
%                               that particular path.
%                               If no type is specified, then the path is
%                               saved to a default type, 'last'. Thus the
%                               only condition is that you can't use 'last'
%                               as your own type!
%                               N.B. Type strings are case-insensitive, so
%                               'ImFiles' is the same as 'imfiles'
%
% Outputs:
%       none
%
% Examples:
%
% Example 1:      
%       % Get the currently defined default directory
%       default_path = retrieve_path()
% 
%       % Navigate to a different directory, and select an arbitrary file
%       [name path] = uigetfile('','Toms_example',default_path)
% 
%       % Update the database for the most recent path
%       update_default_path(path)
% 
%       % Check
%       % Now, run this same example a second time. You'll see that now, 
%       % uigetfile opens up in the same directory as the file you selected
%       % the first time.
%
% Example 2:
%       % Utilises types. In this case, uigetfile will open in the same
%       directory as the most recently accessed file of the specified type,
%       rather than simply the most recently accessed file.
% 
%       % Get default paths for two different types
%       img_path = retrieve_path('imagefiles')
%       snd_path = retrieve_path('soundfiles')
% 
%       % Then select files as necessary
%       [img_name img_path] = uigetfile('','Select an Image File',img_path)
%       [snd_name snd_path] = uigetfile('','Select a Sound File',snd_path)
%       
%       % Update the default path for the types type
%       update_default_path(img_path,'imagefiles')
%       update_default_path(snd_path,'soundfiles')
%
%       % Once this example is run more than once (just select arbitrary
%       % files, for the sake of demonstration) you'll see that img_path
%       % and snd_path are not necessarily the same.
%
%       
% Other m-files required:   retrieve_path.m
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
% Last revised:     16 September 2009

% Set the directory in which you'll write the file
rootdir = which('update_default_path');
rootdir = rootdir(1:end-21);

% Read currently set default values from a .mat file
try
    def_dirs = load([rootdir 'default_paths.mat']);
catch  %#ok<CTCH>
    % If file won't load, no problems. It may have been accidentally 
    % deleted, or not created at all just yet. It will be created when we 
    % save at the end of this function. Use empty catch statement here to avoid 
    % an error if that happens.
end

% Check that it is actually a valid path
if (~ischar(path)) || (exist(path,'dir') ~= 7)
    % Warning, not error, as this is a nonessential function, and I don't
    % want your code to crash and burn on my account!
    warning('update_default_path.m: Invalid path will not be added to defaults database.')
    return
end

% Set the type 
if nargin == 1;     type = 'last';
else                type = lower(varargin{1});
end

% Parse into the structure using type as a dynamic field name
def_dirs.(type) = path;

% Save the structure back to the database file
save([rootdir 'default_paths.mat'],'-struct','def_dirs');
