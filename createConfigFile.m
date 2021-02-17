function createConfigFile( directory, folder, fileName, commentString, ...
    time, radius, angle, velocity, pressureBackground )
%CREATECONFIGFILE Create a configuration file based on model parameters

% Compute temporal, spatial, angular, and velocity limits
timeMin = min(time);
timeMax = max(time);
timeDelta = (timeMax - timeMin) / numel(time);
radiusMin = min(radius);
radiusMax = max(radius);
radiusDelta = (radiusMax - radiusMin) / numel(radius);
angleMin = min(angle);
angleMax = max(angle);
angleDelta = (angleMax - angleMin) / numel(angle);
velocityMin = min(velocity);
velocityMax = max(velocity);
velocityDelta = (velocityMax - velocityMin) / numel(velocity);

% Make directory if none exists
if ~isfolder([directory folder])
    mkdir([directory folder])
end

% Generate file path
FILE_PATH = [directory folder '\' fileName '_' ...
    datestr(now, 'yyyy-mm-dd_HHMM') '.txt'];

% Open file for writing
fileID = fopen(FILE_PATH, 'wt');

% Write information to file
fprintf(fileID, [ ...
    'CONFIGURATION FILE\t' datestr(now, 'yyyy-mm-dd\tHH:MM') '\n\n' ...
    'Temporal limits [s]:\t\t\t'     num2str(timeMin, 3) ...
                               ' : ' num2str(timeDelta, 3) ...
                               ' : ' num2str(timeMax, 3) '\n' ...
    'Radial limits [m]:\t\t\t'       num2str(radiusMin, 3) ...
                               ' : ' num2str(radiusDelta, 3) ...
                               ' : ' num2str(radiusMax, 3) '\n' ...
    'Angular limits [deg]:\t\t\t'    num2str(angleMin, 3) ...
                               ' : ' num2str(angleDelta, 3) ...
                               ' : ' num2str(angleMax, 3) '\n' ...
    'Velocity limits [m/s]:\t\t\t'   num2str(velocityMin, 3) ...
                               ' : ' num2str(velocityDelta, 3) ...
                               ' : ' num2str(velocityMax, 3) '\n' ...
    'Background pressue [mbar]:\t\t' num2str(pressureBackground, 3) '\n' ...
    '\n\tComments:' '\n' commentString]);

fclose(fileID);
end

% 'Oxidation [true/false]:' '\t\t\t' mat2str(Oxidize) '\n' ...
% 'Collision cross-sections [m2]:' '\t\t' mat2str(Sigma_Bg,2) '\t' '[Sr-O2, Ti-O2, O-O2, SrO-O2, TiO-O2, TiO2-O2]' '\n' ...
%{
fprintf(fid, ['Number of type of particles,\t\t' num2str(t_p)  '\n']);
fprintf(fid, ['UC volume,\t\t' num2str(UC_vol)  '\n']);
fprintf(fid, ['Mass atom,\t\t' num2str(m_B)  '\n']);
fprintf(fid, ['Mass Bg,\t\t' num2str(m_O)  '\n']);
fprintf(fid, ['Mass particles,\t\t' num2str(m_P)  '\n']);
fprintf(fid, ['Laser Energy,\t\t' num2str(dis_E_x)  '\n']);
fprintf(fid, ['N_points,\t' num2str(dis_NPoints)  '\n']);
fprintf(fid, ['absorption coefficient,\t\t' num2str(A_coeff)  '\n']);
fprintf(fid, ['binding energy crystal,\t\t' num2str(Eb_UC)  '\n']);
fprintf(fid, ['Density target,\t\t' num2str(Density_target)  '\n']);
fprintf(fid, ['Number of atoms,\t' num2str(N_atoom_totaal)  '\n']);
fprintf(fid, ['CosN,\t\t\t' num2str(cos_theta)  '\n']);
fprintf(fid, ['Degs,\t\t\t' num2str(N_rad)  '\n']);
fprintf(fid, ['Velocity distribution,\t' num2str(sqrt((dis_E_x/(A_coeff*N_atoom_totaal/sum(T_sum))-Eb_UC)*(2/(sum(T_sum)*m_P(n_A+1)))))  '\n']);
fprintf(fid, ['Velocity Width,\t\t' num2str(V_delta)  '\n']);
fprintf(fid, ['E_O_upper,\t' num2str(E_O_high_Sr)  '\n']);
fprintf(fid, ['E_O_lower,\t\t' num2str(E_O_low)  '\n']);
fprintf(fid, ['Spot size,\t\t' num2str(Spot_size)  '\n']);
fprintf(fid, ['3D plot dried?,\t\t' num2str(dried)  '\n']);
fprintf(fid, ['MO2 possible,\t\t' num2str(MO2)  '\n']);
fprintf(fid, ['Vminimaal,\t\t' num2str(V_minimaal)  '\n']);
fprintf(fid, ['Volume BG,\t\t' num2str(4 * pi/3 * (((1)*X_delta)^3-((1-1)*X_delta)^3) * (-cosd(Rad_delta)+cosd(Rad_delta-Rad_delta)))  '\n']);
fprintf(fid, ['Target substrate distance,\t\t' num2str(X_TS)  '\n']);
fprintf(fid, ['maximum number of collisions,\t\t' num2str(limbots)  '\n']);
fprintf(fid, ['Minimal amount of particles,\t\t' num2str(N_local_dist_min)  '\n']);
fprintf(fid, ['smooth_fraction,\t\t' num2str(smooth_fraction)  '\n']);
fprintf(fid, ['smooth_value,\t\t' num2str(smooth_value_set)  '\n']);
%}
