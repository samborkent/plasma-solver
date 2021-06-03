function [bgMatrix, binVolume] = fillBGMatrix( bgDensity, radius, radiusDelta, ...
                                               angle, angleDelta, iAngle )
% function [bgMatrix, binVolume] = fillBGMatrix( bgDensity, radius, angle, angleDelta, iAngle, radiusDelta )
%FILLBGMATRIX Fill background matrix with particles

% Number of radial bins
nRadius = numel(radius);

% Pre-allocate bin volume array
binVolume = zeros(1, nRadius);

for iRadius = 1 : nRadius
    % Calculate bin volumes [m^3]    
    if iRadius == nRadius
        binVolume(iRadius) = (deg2rad(angleDelta) / 3) ...
                     * ( (radius(iRadius) + radiusDelta)^3 - radius(iRadius)^3 ) ...
                     * ( cosd(angle(iAngle)) - cosd(angle(iAngle + 1)));
    else
        binVolume(iRadius) = (deg2rad(angleDelta) / 3) ...
                     * ( radius(iRadius + 1)^3 - radius(iRadius)^3 ) ...
                     * ( cosd(angle(iAngle)) - cosd(angle(iAngle + 1)));
    end                         
end

% Insert number of background particles per radial bin into
%   first velocity bin (v = 0)
bgMatrix = round( bgDensity .* binVolume );

end

