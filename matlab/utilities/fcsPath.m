function [path] = fcsPath()
%FCSPATH Returns path of the FCS on the present system

path = fileparts(which('fcsPath'));

% It's in the utilities directory, so lose the trailing /utilities
path = path(1:end-9);

end

