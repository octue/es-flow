function [verStr, gitInfo] = fcsVer(flag)
%FCSVER returns current version (or master-committed datestamp) of the FCS
% Checks that the current version of the FCS is synced with the remote 
% repo, and returns the datestamp of the last change to the remote repo
% along with a branch name. In this way, the branch and datestamp of exactly
% what FCS version is running can always be retrieved. Fails with error if the
% local copy is not contained in the remote origin.
%
% Syntax:
%       verStr = fcsVer Returns a human readable string(the remote repo commit
%       count) corresponding with a particular git SHA1 code in a remote repo
%
%       [verStr, gitInfo] = fcsVer Also returns a structure containing long and
%       short form Hashes of the current working copy.
%
%       [verStr, gitInfo= = fcsVer('-warn') Throws a warning, instead of a hard
%       error, when the local repository is not clean, commited, and pushed to
%       remote.
%
% Inputs:
%       none
%
% Outputs:
%
%       verStr          string      Version of the FCS, 'v_****', where ****
%                                   is the remote repo commit counter for the
%                                   current clean, pushed, working copy. In
%                                   -warn mode when the repo is uncommitted, the
%                                   string is simply 'v_uncommitted'.
%
%       gitInfo         struct      Structure with the following fields:
%                   .ShortStatus        string containing the result of git status -s
%                   .Hash               string containing the git hash of the
%                                       current working copy
%                   .ShortHash          string the short has equivalent
%                   .Count              integer commit count of the remote commit
%                                       corresponding to the current 
%                                       working copy hash. This gives, in
%                                       effect, a human-readable version number.
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
% Revision History:        	13 January 2016     Created from tssVer()
%
% Copyright (c) 2016 Ocean Array Systems, All Rights Reserved.

% Store the current working directory for later
workDir = pwd;

% Change to the current FCS directory
cd(fcsPath)

try
    % Analyse git properties and store in structure
    [success, gitInfo.ShortStatus] = system('git status -s');
    assert(success==0,'FAILED: git status -s')
    assert(isempty(gitInfo.ShortStatus), 'Local copy not clean. Commit all changes and push to remote.')

    [success, gitInfo.Hash] = system('git rev-parse HEAD');
    assert(success==0,'FAILED: git rev-parse HEAD')
    [success, gitInfo.ShortHash] = system('git rev-parse --short HEAD');
    assert(success==0,'FAILED: git rev-parse --short HEAD')
    [success, count] = system('git rev-list HEAD --count');
    gitInfo.Count = str2double(count);
    assert(success==0,'FAILED: git rev-list HEAD --count')

    % Check if the present commit is in the remote repo
    [out, txt] = system(['git branch -r --contains ' gitInfo.Hash ' origin/*']);
    assert((out == 0) && ~isempty(txt), 'Present local copy hash not available in remote repo. Commit all changes and push to remote.')
        
    % Set the useable version string as the commit count
    verStr = ['v_' num2str(gitInfo.Count)];
    
catch me
    
    % Issue a warning instead of a hard error if there's an input '-warn' flag
    if nargin >= 1 
        if strcmpi(flag,'-warn')
	         warning('Local copy not clean, or present local copy hash not available in remote repo. Commit all changes and push to remote.')
        end
        verStr = 'v_uncommitted';
        cd(workDir)
    else
        throw(me)
    end
    
end

% Change back to working directory
cd(workDir)
