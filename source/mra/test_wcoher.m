% Load the data
testMARS


% Establish the scale range from 0 to 50 (is this dyadic scale? if not then
% we need a way wider range)
scales = 1:512;

ntw = 4; % Number of time windows used for smoothing CWT coefficients before computing WCOH and WCS
wname  = 'cmor1-3';

wcoher(ux(:),uy(:),scales,wname,'ntw',4,'plot','cwt');