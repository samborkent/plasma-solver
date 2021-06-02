function [bgMatrix, binVolume] = fillBGMatrix( bgDensity, radius, angleDelta )
% function [bgMatrix, binVolume] = fillBGMatrix( bgDensity, radius, angle, angleDelta, iAngle, radiusDelta )
%FILLBGMATRIX Fill background matrix with particles

nRadius = numel(radius);

binVolume = zeros(1, nRadius);

for iRadius = 1 : nRadius-1
    % Calculate bin volumes [m^3]    
%     if iRadius == nRadius
%         binVolume(iRadius) = (4/3)*pi ...
%                                 * ( (radius(nRadius) + radiusDelta )^3       ...
%                                 -   radius(iRadius)^3 ) ...
%                                 * ( cosd(angle(iAngle))     ...
%                                 -   cosd(angle(iAngle + 1)) );
%     else
%         binVolume(iRadius) = (4/3)*pi ...
%                                 * ( radius(iRadius + 1)^3       ...
%                                 -   radius(iRadius)^3 ) ...
%                                 * ( cosd(angle(iAngle))     ...
%                                 -   cosd(angle(iAngle + 1)) );
%     end
    binVolume(iRadius) = (1/3)*deg2rad(angleDelta)*(1 - cosd(angleDelta)) ...
                         * ( radius(iRadius + 1)^3 - radius(iRadius)^3 );
end

% Insert number of background particles per radial bin into
%   first velocity bin (v = 0)
bgMatrix = bgDensity .* binVolume;

end

