function disperror(err)
%DISPLAYERROR Displays a formatted error message to the console
%   DISPLAYERROR takes an error structure as an input and displays the 
%   error to the console in a manner similar to MATLAB.  The output 
%   provides links to the files containing the error as well as the line
%   text where the error occurred.
%
%   DISPLAYERROR displays the lasterror if no input is specified.
%
%    Examples:
%  
%     displayError(err)     Displays the error described by the error
%                           structure err
%
%     displayError          Displays the last error as described by
%                           lasterror
%  
%   See also lasterror, error, try, catch, rethrow.
%
%   Reference page in Help browser
%      <a href="matlab:doc displayError">doc displayError</a>

%   Written by Chris J Cannell
%   Contact ccannell@gmail.com for questions or comments.
%   02/13/2007

% input argument error checking
if nargin == 0
    err = lasterror;
elseif nargin > 1
    error('Number of input arguments must be zero or one');
end

% check that input is an error structure
if ~isa(err,'MException') || ~isprop(err, 'message') || ~isprop(err, 'identifier') || ~isprop(err, 'stack')
    error('Input is not a valid error structure');
end

msgString = getReport(err,'extended','hyperlinks','on');
disp(msgString)
% % display error message
% fprintf('\n??? %s\n\n', err.message);
% 
% % display list of files and line numbers were error occurred
% for i = 1:length(err.stack)
%     fprintf('Error in ==> <a href="error:%s,%d,1">%s at %d</a>\n',...
%         err.stack(i,1).file, err.stack(i,1).line, err.stack(i,1).name,...
%         err.stack(i,1).line);
%     fileLine = getLineFromFile(err.stack(i,1).file,err.stack(i,1).line);
%     fprintf('%s\n\n', fileLine);
% end
