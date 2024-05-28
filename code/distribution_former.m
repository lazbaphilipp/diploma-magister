function [distributionValues,phases] = distribution_former(dx,f, phi,amp_dB)
if(size(phi) ~= size(amp_dB))
    error("sizes unequal")
end

amp = db2pow(amp_dB);
distributionValues = zeros(size(dx));
phases = zeros(length(dx),length(phi));
for ind = 1:length(phi)
phases(:,ind) = 2.*pi.*f.*dx.*sind(phi(ind))/physconst("Lightspeed");
distributionValues = distributionValues + amp(ind).*exp(1i.*phases(:,ind));
end

end