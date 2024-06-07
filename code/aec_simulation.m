function [channel_res,theta] = aec_simulation(dVals,dx,f, resolution, theta)
if nargin < 4
    resolution = 0.5;
end
if nargin < 5
    theta = (-90:resolution:90);
end
channel_coeff = exp(1i.*(-2.*pi.*f.*dx.*sind(theta)/physconst("Lightspeed")));
channel_res = sum(dVals.*channel_coeff,1);
% theta = theta';
end