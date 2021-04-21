clear; close all; clc;

%--------------------------------------------------------------------------
% Dimensional limits
%--------------------------------------------------------------------------

% Angular limits
angleMin    = 0;        % Start angle [deg]
angleMax    = 90;       % End angle [deg]
angleDelta  = 3;        % Radial step size [deg]

% Temporal limits
timeMin     = 0;        % Start time [s]
timeMax     = 10E-6;    % End time [s]
timeDelta   = 1E-7;     % Time step duration [s]

% Radial limits
radiusMin   = 0;        % Start position [m]
radiusMax   = 0.06;     % End position [m]
radiusDelta = 10E-5;     % Radial step size [m]

% Velocity limits
veloMax         = 2.5E4; % Maximal initial velocity
veloStepsize    = 1;     % Velocity step size

%--------------------------------------------------------------------------
% Axis creation
%--------------------------------------------------------------------------

% Angular axis
angle   = angleMin : angleDelta : angleMax - angleDelta;
nAngle  = numel(angle);

% Temporal axis
time    = timeMin : timeDelta : timeMax - timeDelta;
nTime   = numel(time);

% Radial axis
radius  = radiusMin : radiusDelta : radiusMax - radiusDelta;
nRadius = numel(radius);

% Velocity array
veloDelta   = (radiusDelta / timeDelta) / veloStepsize;
velo        = 0 : veloDelta : veloMax - veloDelta;
nVelo       = numel(velo);

%--------------------------------------------------------------------------
% Additional
%--------------------------------------------------------------------------

% MO: 4, ABO: 8, ABCO: 11

nSpecies = 11;

% N: 1, k: 2

nField = 2;

%--------------------------------------------------------------------------
% Matrix creation
%--------------------------------------------------------------------------

nLoops = 0;

matrix = floor( randn(nVelo, nSpecies, nRadius) + 1.1 );

matrixBG = floor( randn(nVelo, nRadius) + 1.1 );

tic;

for iAngle = 1 : 1 % nAngle
    
for iTime = 1 : nTime

% Only loop through radial bins where particles could potentially be
%   present
if (nVelo * iTime) < nRadius
    endRadius = (nVelo * iTime) - 1;
else
    endRadius = nRadius - 1;
end
    
for iRadius = endRadius : -1 : 1
    
for iVelo = nVelo : -1 : 2
        
    A = tic;
    
    nPlasma = matrix(iVelo, :, iRadius);
    
    B = toc(A);

%     % Skip bins
%     if nPlasma < 1
%         continue
%     end
% 
%     nTemp = nPlasma;

%     % Restrict number of traveled bins to the total number of
%     %   radial bins
%     if ( iRadius + iVelo - 1 ) > nRadius
%         nRadiusDelta = nRadius - iRadius;
%     else
%         nRadiusDelta = iVelo - 1;
%     end
% 
%     for jRadius = 0 : (nRadiusDelta - 1)
% 
%     if nTemp < 0
%         break
%     end
% 
%     thisRadius = iRadius + jRadius;
% 
%     for jVelo = 1 : (iVelo - 1)
% 
%         if nTemp < 0
%             break
%         end
% 
%         nBG = matrixBG(jVelo, iRadius + jRadius);
% 
%         if nBG < 1
%             continue
%         end
% 
%         nCol = nBG + nPlasma;
% 
%         if nCol < 3
%             continue
%         end
% 
%         nTemp = nTemp - 0.001;
% 
%         nLoops = nLoops + 1;
% 
%     end % jVelo
% 
%     end % jRadius
    
end % iVelo

end % iRadius
    
end % iTime
    
end % iAngle

toc;

% x = myfun(x);
% 
% function x = myfun(x)
% x = 1.2*x;

% N = 3e3;
% 
% x = randn(N);
% x = x * 1.2;

% N = 2e3;
% 
% x = randn(N);
% y = zeros(N);
% 
% for col = 1:N
%     for row = 1:N
%         if x(row, col) > 0
%             y(row, col) = x(row, col);
%         end
%     end
% end

% N = 10e3;
% 
% x = zeros(N, 1);
% 
% x(1) = 1000;
% 
% for k = 2:N
%     x(k) = 1.05*x(k-1);
% end