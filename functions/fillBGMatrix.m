function [bgMatrix, binVolume] = fillBGMatrix( bgDensity, nVelo, radius, angle, iAngle )
%FILLBGMATRIX Fill background matrix with particles

nRadius = numel(radius);

binVolume = zeros(1, nRadius);

bgMatrix = zeros(nVelo, nRadius);

for iRadius = 1 : nRadius
    % Calculate bin volumes [m^3]
    if iRadius == 1
        binVolume(iRadius) = (4/3)*pi ...
                                * radius(iRadius)^3     ...
                                * ( cosd(angle(iAngle)) ...
                                -   cosd(angle(iAngle + 1)) );
    else
        binVolume(iRadius) = (4/3)*pi ...
                                * ( radius(iRadius)^3       ...
                                -   radius(iRadius - 1)^3 ) ...
                                * ( cosd(angle(iAngle))     ...
                                -   cosd(angle(iAngle + 1)) );
    end

    % Insert number of background particles per radial bin into
    %   first velocity bin (v = 0)
    bgMatrix(1, iRadius) = bgDensity * binVolume(iRadius);
end

end

