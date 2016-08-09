function [power] = herbieYieldIntegrator(params, userData)
%HERBIEYIELDINTEGRATOR Takes a power curve in the user data, and velocity as the
%input parameter. It returns the power generated. When sampled this way in
%herbie using a velocity distribution PDF, this is the average power across time
%so simply multiply by the number of seconds in the time period to get total
%yield.
% Note that the power curve reference velocity must be the same reference
% velocity as the distribution.
%
% TODO create a power curve object encompassing this use case of herbie as a
% method, and managing the reference base.

if isa(userData,'PowerCurve')
    % TODO implement power curve lookup method
else
    % Userdata is a structure containing velocity and power fields
    power = interp1(userData.velocity, userData.power, params(1), 'linear');
    if isnan(power)
        disp(params)
        
        error('nan')
    end
end

