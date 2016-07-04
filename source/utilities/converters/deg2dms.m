function [out] = deg2dms(in)
%DEG2DMS Convert lists of values in decimal degrees to values in degrees,
% minutes, seconds

Lat = in(:,1);
Lon = in(:,2);

if size(Lat,2)>1 || size(Lon,2)>1
    error('Lat and Lon must be column vectors')
end

% Round toward zero
LatDegrees = fix(Lat);
LonDegrees = fix(Lon);

% Decimal minutes (absolute value)
LatDM = 60*abs(Lat-LatDegrees);
LonDM = 60*abs(Lon-LonDegrees);

% Round toward zero
LatMinutes = fix(LatDM);
LonMinutes = fix(LonDM);

% Decimal seconds (absolute value)
LatDS = 60*abs(LatDM-LatMinutes);
LonDS = 60*abs(LonDM-LonMinutes);

% Compile into matrices
LatDMS = [LatDegrees LatMinutes LatDS];
LonDMS = [LonDegrees LonMinutes LonDS];

out = [LatDMS LonDMS];

fprintf('%11.0f, %11.0f, %11.6f, %11.0f, %11.0f, %11.6f\n',out')

end

