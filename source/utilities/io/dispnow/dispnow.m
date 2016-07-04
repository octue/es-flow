function dispnow(string)
%DISPNOW Displays output to the command line with a datestamp, as per the
%requirement for TurbineGrid outputs.
disp([datestr(now) '    ' string])