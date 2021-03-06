%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plasma Solver rad version v3.1.5 03-05-20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Made by Tom Wijnands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sub scripts
%   -Makestruct_v2              - generate the initial situation
%   -Plot_3d_full            - Reshape space from 1D to 3D
%   -Plot_Arr_full           - Calculations the arriving species
%   -myfastint2str           - fast number to string conversions
%
%   ---General layout---
%   input settings
%   calculations per angle (for loop)
%       generate space matrix & initial velocities
%       calculations per timestep (for loop)
%           caculations per distance step/container x-axis (for loop)
%               collision probability calculations + calculate moved distance
%               check if particles arrive at substrate
%               check and oxidize particles if needed
%               particles moved to new location
%           smooth and normlize        
%           plot 1d graphs     
%       plot 1d angle plots
%   calculate and plot 3d images
%   sum and plot all arriving particles
%

                                
% typical naming conventions 
%"X"        -for distance calcuations
%"N"        - for partilces calcuations
%"V"        - for veolicty calcuations
%"K"        - number of collisions
%"P"        - type of particle
%"mo"       = moving particles calcuations
%"st"       - standing partilces ( not moving) calcuations
%"bg"    - are for the calculations of the background.
%"bots"     - colliding particles
%"nbots"    - not colliding particles 
%"M"        - current space matrix of all plume particles
%"M_1"      - new  space matrix of all plume particles
%"Mbg"      - current space matrix of all background particles
%
% examples:
% x_new_bots_mo, distance/new/colliding particles/ with moving particles, the new distance of colliding particles with a moving background
% x_new_bots_st, distance/new/colliding particles/with standing particles,  the new distance of colliding particles with a standing background
%
% in the output folder structure:
% t_p_0 -> all particles
% t_p_1 -> M particle as in Ti in TiO2
% t_p_2 -> O particle as in O in TiO2
% t_p_3 -> MO particle as in TiO in TiO2
% t_p_4 -> MO2 particle as in TiO2 in TiO2


 


clear all
clearvars -except files nn
close all
profile off % on  for optimalization testing
warning('OFF')


%% init starting conditions
%-------------------------------------------------------------------------
% Simulation settings
%---------------------------------------------------------

%Save location
save_loc = 'C:\Users\Sam\Google Drive\School\Master\Capita selecta\Model_Tom\TiO2_test_data\'; %top domain
folder = 'TiO_2_simulation_xx_mBar\'; % This measurement folder
plek = [save_loc folder]              % location

comment=('F_xxTiO2');  %make a comment on the simulation, saved in Zdata.txt

% plots
dried = 0;                          %make 3d images 1=yes 0=no
tijden = [5,15,20,30];              %the time values to make 3d plot
Arr_sum = 0;                        %Calculate arriving particles  1=yes 0=no
make_1d_plots = 1 ;                % make 1d plots at start rad

%limitations
MO2 = 1;                            %Let particles form MO2,  1=yes 0=no
N_local_dist_min = 1;               %minimal amount of particles
limbots = 7;                        %maximum number of collisions possible
Min_angle_particles = 1E7;          %Do not start calculations at angles with less than these number of particles
show_BG_inplot = 0;                 %show background gas in plot, (normal off, 0 )
smooth_fraction = 1E4;              %fraction to smooth
smooth_value_set = 15;                  %standaard matlab moving averige smooth value

%---------------------------------------------------------

%
% Simulation enviroment settings
%---------------------------------------------------------
%time settings
T_begin = 0;                %Starting time of the simulation in [s]
T_end = 18E-6;              %End time of the simulation in [s]
T_delta = 1E-7;             %Timestep in [s]

% distance settings
X_begin = 0;                %begin distance in [m]
X_end = 0.06;               %end X distance in [m]
% X_delta = 0.5E-4;           %Distance stepsize in [m]
X_delta = 0.005;

% angle settings
Rad_start = 0;              %Start angle in [deg]
Rad_stop = 90;
N_rad = 90;                 %number of total angles in [deg]
% Rad_delta = 3;              %angle_step  in [deg]
Rad_delta = 30;


% Deposition settings
P_dep = 0.1;                %depostion pressure in [mBar]
X_TS = 0.05;                %target substrate distance in [m]
dis_E_x = 0.052;            %laser energy in [J]
Spot_size = 2*10^(-3)*1*10^(-3)*100*10^(-9); %spot size Length*Width*Depth in [m2]
kb = 1.38E-23;              %boltzman constant
Temperature_BG = 300;        %Temperature background in [K]
density_bg = P_dep*100/kb/Temperature_BG;    %density bacground, ideal gas law

%---------------------------------------------------------

% Material settings 
%---------------------------------------------------------
%AOx
UC_vol = 4.5937E-10^2*2.9587E-10; %volume AOx UC a*b*c
Density_target = 1;               %Density_target
A_coeff = 0.66;                   %adsorbtion coeff AOx 1=100% absorbed
m_c = 7.948E-26;                  %mass A atom
m_o = 1.32843216E-26*4;           %mass O2 molecule
Eb_UC = 3.254*1.6E-19;            %binding energy crystal in J

Oxidize = 1 ;                    %let the particle oxidize, 1 = no oxidization possible
E_O_high = 200*1.602176634*10^(-19);%3.2E-17*Oxidize; %200*1.6E-19;      %upper bound energy oxidation %als 0 geen oxidatie
E_O_low = 0.5*1.6E-19;     %lower bound energy oxidation
T_sum = [1 2];                    %amount of atoms per unit cell
t_p = 3;                          %number of start particles +1
t_o = 2;                          %number of oxydation states
t_oo = 2;                         %het oxide deeltje plaats.
t_pp = 5;                         %total number of particles
m_p = [m_c m_o/2 m_c+m_o/2 m_c+m_o]; %A O AO AO2
%cross section background gas molecule
Opp_plasma=[pi*(17.6E-11+4.8E-11*2)^2 pi*(3*4.8E-11)^2 pi*(17.6E-11+4.8E-11*3)^2 pi*(17.6E-11+4.8E-11*4)^2];   % [ Ti-O2, O2-O2, TiO-O2, TiO2-O2 ]

%---------------------------------------------------------

% Initial plume settings
%---------------------------------------------------------
V_delta = 1500;                     %Width velocity distribution
cos_theta = 25;                     %radial sharpe of the plasma
dis_NPoints =  100;                 %Number of datapoints in the initial velocity distribution

%---------------------------------------------------------
% % Setting up the initial distributions
%---------------------------------------------------------
N_atoom_totaal = Spot_size/UC_vol*Density_target/sum(T_sum);    %Amount of ablated atoms per UC
N_atoom_totaal_N = T_sum*N_atoom_totaal;                          %total amount of ablated atoms
theta = 0:Rad_delta:N_rad;                                      % the angles

for theta_i = 1:size(theta')-1
    N_atom_temp(theta_i) = ((4*pi/3 .* (X_delta^3) .* (-cosd(theta(theta_i+1))+cosd(theta(theta_i)))).*(cosd(theta(theta_i))).^cos_theta);
end
N_atoom_rad = N_atom_temp./sum(N_atom_temp).*N_atoom_totaal;
V_minimaal = X_delta/T_delta;                 %Minimal velocity to move forward one container
V_maximaal = 40000;                           %maximal initial velocity
V_stepsize = 5;                               %Initiele snelheid dv 0:(V_minimaal/V_stepsize):V_maximaal;

%---------------------------------------------------------
%write settings into data config file
%---------------------------------------------------------
mkdir([save_loc folder])
fid = fopen([save_loc folder '\ZData.txt'],'wt');
fprintf(fid, ['T_begin,\t\t' num2str(T_begin)  '\n']);
fprintf(fid, ['T_end,\t\t\t ' num2str(T_end)  '\n']);
fprintf(fid, ['T_delta,\t\t' num2str(T_delta)  '\n']);
fprintf(fid, ['X_begin,\t\t' num2str(X_begin)  '\n']);
fprintf(fid, ['X_end,\t\t\t' num2str(X_end)  '\n']);
fprintf(fid, ['X_delta,\t\t' num2str(X_delta)  '\n']);
fprintf(fid, ['Rad_delta,\t\t' num2str(Rad_delta)  '\n']);
fprintf(fid, ['P_dep,\t\t' num2str(P_dep)  '\n']);
fprintf(fid, ['Pressure [mBar],\t\t' num2str(density_bg/100*kb*300)  '\n']);
fprintf(fid, ['Oxidize possible?,\t\t' num2str(Oxidize)  '\n']);
fprintf(fid, ['Cross section,\t\t' num2str(Opp_plasma)  '\n']);
fprintf(fid, ['Number of type of particles,\t\t' num2str(t_p)  '\n']);
fprintf(fid, ['UC volume,\t\t' num2str(UC_vol)  '\n']);
fprintf(fid, ['Mass atom,\t\t' num2str(m_c)  '\n']);
fprintf(fid, ['Mass Bg,\t\t' num2str(m_o)  '\n']);
fprintf(fid, ['Mass particles,\t\t' num2str(m_p)  '\n']);
fprintf(fid, ['Laser Energy,\t\t' num2str(dis_E_x)  '\n']);
fprintf(fid, ['N_points,\t' num2str(dis_NPoints)  '\n']);
fprintf(fid, ['absorption coefficient,\t\t' num2str(A_coeff)  '\n']);
fprintf(fid, ['binding energy crystal,\t\t' num2str(Eb_UC)  '\n']);
fprintf(fid, ['Density target,\t\t' num2str(Density_target)  '\n']);
fprintf(fid, ['Number of atoms,\t' num2str(N_atoom_totaal)  '\n']);
fprintf(fid, ['CosN,\t\t\t' num2str(cos_theta)  '\n']);
fprintf(fid, ['Degs,\t\t\t' num2str(N_rad)  '\n']);
fprintf(fid, ['Velocity distribution,\t' num2str(sqrt(((dis_E_x*0.9*A_coeff)/3/N_atoom_totaal-Eb_UC)/2/m_p(1)))  '\n']);
fprintf(fid, ['Velocity Width,\t\t' num2str(V_delta)  '\n']);
fprintf(fid, ['E_O_upper,\t' num2str(E_O_high)  '\n']);
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
fprintf(fid, ['Comments,\t\t' comment  '\n']);
fclose(fid);



%% Begin calculations for all angles
%--------------------------------------------------------------------------

for rad = Rad_start+Rad_delta: Rad_delta : Rad_stop
    
    %% Setup initial conditions
    %--------------------------------------------------------------------------
    
    %pre allocate space for saving time
    Totaal_t = zeros(round((X_end-X_begin)/X_delta),round((T_end-T_begin) / T_delta));
    N_max_2 = zeros(round((T_end-T_begin) / T_delta),limbots+1);
    N_loc_2 = zeros(round((T_end-T_begin) / T_delta),limbots+1);
    X_str_end1 = ['X_' myfastint2str(((X_end-X_begin)/X_delta)-1)];
    X_str_end2 = ['X_' myfastint2str(((X_end-X_begin)/X_delta)-2)];
    rad
    rad_string = myfastint2str(rad);
    N_atoom = N_atoom_rad(rad/Rad_delta);         %Number of atoms at this angle
    % Skip angles with lower then Min_angle_particles particles
    if N_atoom < Min_angle_particles
        N_atoom = 1;
    end
    %create folder to save the images of the results
     make_angle_folder_structure_Tom(save_loc, folder, t_p, t_o,rad) % make folder structure
     
    %% Calculate the initial velocity praticle distribution
    %--------------------------------------------------------------------------
    % set average and width of the gaussian distribution
    % make the inital N-V distribution
    for mm = 1:(t_p-1)
        clear dis_v_dist_x dis_N_dist_x
        v_gem = sqrt(((dis_E_x*0.9*A_coeff)/3/N_atoom_totaal-Eb_UC)/2/m_p(mm));
        dis_vPlotMax = v_gem*2.5;     % Maximum speed to plot
        dis_v_dist_x = 0:(V_minimaal/V_stepsize):V_maximaal;
        for dis_iPoint = 1:length(dis_v_dist_x)
            %dis_v_dist_x(dis_iPoint) = (dis_iPoint-1)/(dis_NPoints) * (dis_vPlotMax);
            dis_N_dist_x(dis_iPoint) = 1/(V_delta*sqrt(2*pi)) * exp(-(dis_v_dist_x(dis_iPoint)-v_gem)^2/(2*V_delta^2));
        end
      
        nor = ((N_atoom*T_sum(mm))/sum(sum(dis_N_dist_x))); %normalizatie the number of atoms
        %plot of the initial velocity distribution and save
        if rad == Rad_delta
            fig = figure('visible','on');
            plot(dis_v_dist_x,dis_N_dist_x.*nor,'x--')
            xlabel('Velocity [m/s]')
            ylabel('Number of particles')
            set(gca,'FontSize',16)
            set(gca,'LineWidth',2)
            set(get(gca,'xlabel'),'FontSize',18)
            set(get(gca,'ylabel'),'FontSize',18)
            set(get(gca,'title'),'FontSize',18)
            set(findobj('Type','line'),'MarkerSize',10)
            set(findobj('Type','line'),'LineWidth',3)
            if rad < 10
                saveas(fig,[save_loc folder '\t_p_' myfastint2str(mm) '\Rad0' rad_string '\img\Init_Velocity_' myfastint2str(mm)   '.jpg'])
                saveas(fig,[save_loc folder '\t_p_' myfastint2str(mm) '\Rad0' rad_string '\img\Init_Velocity_' myfastint2str(mm)   '.fig'])
            else
                saveas(fig,[save_loc folder '\t_p_' myfastint2str(mm) '\Rad' rad_string '\img\Init_Velocity.jpg'])
                saveas(fig,[save_loc folder '\t_p_' myfastint2str(mm) '\Rad' rad_string '\img\Init_Velocity.fig'])
            end %end of  rad < 10
%             close
        end %end if t==1
        %sum(dis_N_dist_x)*nor
        init_N_dist((mm-1)*dis_iPoint+1:mm*dis_iPoint) = dis_N_dist_x.*nor;
        init_V_dist((mm-1)*dis_iPoint+1:mm*dis_iPoint) = dis_v_dist_x;%
        p((mm-1)*dis_iPoint+1:mm*dis_iPoint)=mm-1;
    end %mm=1:(t_p+t_o-1)-1
        
    %% Setup make matrix for calulations
    %--------------------------------------------------------------------------
    % Plasma matrix M and background matrix Mbg
    % plasma matrix t+1 is M_1 and background t+1 matrix Mbg
    k = zeros(1,length(init_V_dist)); %set all collisions k to zero
    M = Makestruct_v2( X_begin , X_end ,  X_delta , 1 , rad , Rad_delta ); %make the matrix in this function
    Total_Ek_theory = sum(m_c/2*dis_v_dist_x.^2.*dis_N_dist_x.*nor); %calcultate the initial sum of kinetic energy
    Total_p_theory = sum(m_c*dis_v_dist_x.*dis_N_dist_x.*nor) ; %calculate the initial momentum
    %initiele condities geven voor de plasma in de matrix
    M = setfield(M,'X_1','V',init_V_dist);
    M = setfield(M,'X_1','N',init_N_dist);
    M = setfield(M,'X_1','K',k);
    M = setfield(M,'X_1','P',p);
    clear x y
    %Calculate the background matrix
    Mbg = Makestruct_v2( X_begin , X_end ,  X_delta , 0 , rad , Rad_delta );
    %write values into matrix
    for x = 1:round((X_end-X_begin)/X_delta);
        Vol_xy = 4 * pi/3 * (((x)*X_delta)^3-((x-1)*X_delta)^3) * (-cosd(rad)+cosd(rad-Rad_delta));
        Nbg_xy(x) = density_bg*Vol_xy;
        Vol_T(x) = Vol_xy;
        %schrijf N in matrix Mbg
        X_n = ['X_' myfastint2str(x)];
        Mbg = setfield(Mbg,X_n,'N',[Nbg_xy(x) 0]);
        Mbg = setfield(Mbg,X_n,'VOL',Vol_xy);
    end
    clear x y X_n k Vol_xy Nbg_xy
    %sum up all the background molecules
    for x = 1:floor((X_end-X_begin)/X_delta)
        X_n = ['X_' myfastint2str(x)];
        O2_bg(x) = sum(Mbg.(X_n).N);
    end
    Total_O2 = sum(O2_bg); %totaal numbers of moleculen in the background
    %% main programma per timestep
    %-------------------------------------------------------------------------
    for t=1:(T_end-T_begin) / T_delta %for all timesteps
        %pre allocate space
        Ek_arr_N = zeros(round((X_end-X_begin)/X_delta),t_p+1);                          %aantal,type
        tijd_string = myfastint2str(t);
        
        %initialize the (t+1) structure
        M_1 = Makestruct_v2( X_begin, X_end,  X_delta, 1, rad, Rad_delta);
        x_names = fieldnames(M); %number of X step in the strucute
        
        for xx = 1:size(x_names) %for loop per X
            %calculate simulation from back to forward in distance
            x = length(x_names)-xx+1;
            ssx = length(fieldnames(M.(x_names{x})));
            if  ssx > 1 %check if empty
                %Get values from the struct M and M_bg for calulations
                %Get information from the structure for the plasma
                V = (M.(x_names{x}).V)';%waardes V  [velocity]
                X = (M.(x_names{x}).X); %waardes X  [Distance]
                N = (M.(x_names{x}).N); %waardes N  [Number of particles]
                K = (M.(x_names{x}).K); %waardes K  [Number of collisions]
                P = (M.(x_names{x}).P); %waardes P  [Type of particle]
                %Get information from the structure for the background
                VOLbg = (Mbg.(x_names{x}).VOL); %waardes VOL
                Kbg = (Mbg.(x_names{x}).K); %waardes K
                Vbg = (Mbg.(x_names{x}).V)';%waardes V
                Nbg = (Mbg.(x_names{x}).N); %waardes N
                size_V_dist = length(V);    %The number of different groups in this container.
                [V_sort_val V_sort_loc] = sort(V,'descend'); %sort on fasted velocity first
                %% Caculations on the individial NV groups
                %-------------------------------------------------------------------------
                for vv = 1:size_V_dist
                    %select the highest velocity first
                    
                    jj = V_sort_loc(vv);
                    %remove isnan from distribution, Should not be necessary
                    if  isnan(V(jj)) %
                        % if the goup has particles transfer them.
                        if round(V(jj)/(X_delta/T_delta))*(X_delta/T_delta)==0 && M.(x_names{x}).N(jj) > N_local_dist_min
                            M_1.(x_names{x}).V(jj) = M.(x_names{x}).V(jj);
                            M_1.(x_names{x}).N(jj) = M.(x_names{x}).N(jj);
                            M_1.(x_names{x}).K(jj) = M.(x_names{x}).K(jj);
                        else
                            N_local_dist=0;
                        end
                    else % if the group has a verlocity
                        N_local_dist = N(jj); %Local distribution of particles
                        V_local_dist = V(jj); %local velocity
                        K_local = K(jj);      % Local number of collisions
                        P_local = P(jj);
                        V_local_dist = round(V_local_dist/(X_delta/T_delta))*(X_delta/T_delta); %rond to neaserst step size.
                        Vbg = round(Vbg/(X_delta/T_delta))*(X_delta/T_delta); %rond to neaserst step size.
                        
                        
                        %Calculte the average density
                        %-----------------------------------------------------------
                        %set some intial condictions before calculating per V
                        
                        S_exp = 0;
                        N_temp_left = 0;
                        R_kans_botsen_st_exp = 0;
                        R_kans_botsen_mo_exp = 0;
                        R_kans_bots_st = 0;
                        R_kans_bots_mo = 0;
                        x_newbg_bots_st = 0;
                        Vbg_exp = 0;
                        Vbg_min_mo= 0 ;
                        
                        S_exp = abs(V_local_dist*T_delta/X_delta);  %Expected distance traveled
                        N_temp_left = N_local_dist;
                        
                        %check a particle travles a distance
                        if S_exp == 0 || V_local_dist == 0 || N_local_dist < N_local_dist_min  || S_exp == Inf
                            % if no distance traveld
                            % set values to 0
                            R_kans_botsen_st_exp = 0;
                            R_kans_botsen_mo_exp = 0;
                            R_kans_bots_st = 0;
                            R_kans_bots_mo = 0;
                            Vbg_min_mo = 0;
                            x_new_bots_st = 0;
                            x_new_bots_mo = 0;
                            x_newbg_bots_st = 0;
                            x_newbg_bots_mo = 0;
                            S_exp = 1;
                            Vbg_min_st = 0;
                            Sbg_st = 0;
                        else
                            %if distance is expected to be traveled
                            for n_nbg = 1:S_exp
                                %calculate per X_delta the colision probability
                                % if the distance is smaller then the step size
                                if x+n_nbg-1 > ((X_end-X_begin)/X_delta)
                                    R_kans_botsen_st_exp(n_nbg) = 0;
                                    R_kans_botsen_mo_exp(n_nbg) = 0;
                                    N_temp_left(n_nbg+1) = N_temp_left(n_nbg);
                                    R_kans_bots_st(n_nbg) = 0;
                                    R_kans_bots_mo(n_nbg) = 0;
                                    Vbg_min_mo(n_nbg) = 0;
                                    S_st(n_nbg) = 0;
                                    S_mo(n_nbg) = 0;
                                    Sbg_st(n_nbg) = 0;
                                    Sbg_mo(n_nbg) = 0;
                                else
                                    % Left off
                                    try
                                        N_loc = M.(x_names{x+n_nbg-1}).N;
                                    catch
                                        N_loc(n_nbg) = 0;
                                    end
                                    %calculatug amount of particles in background
                                    Nbg_exp_V = (Mbg.(x_names{x+n_nbg-1}).N)/(Mbg.(x_names{x+n_nbg-1}).VOL);%denisty background
                                    Nbg_exp_N = (Mbg.(x_names{x+n_nbg-1}).N);%number N atoms
                                    %calculating ht eamount of partciels in the plasma
                                    Nbg_exp = Nbg_exp_N; %+sum(Np_exp);
                                    %calculate colision probability
                                    R_kans_botsen_st_exp(n_nbg) = Nbg_exp_V(1).*N_temp_left(n_nbg)/(sum((N_loc))+N_temp_left(n_nbg))*(Opp_plasma(P_local+1))*X_delta; %only valid if density plasma is lower than background?
                                    R_kans_botsen_mo_exp(n_nbg) = Nbg_exp_V(2).*N_temp_left(n_nbg)/(sum((N_loc))+N_temp_left(n_nbg))*(Opp_plasma(P_local+1))*X_delta;
                                    %R_kans_botsen_mo_exp(n_nbg)=0;
                                    nanis = find(isnan(R_kans_botsen_st_exp))';
                                    R_kans_botsen_st_exp(nanis) = 0;
                                    nanis = find(isnan(R_kans_botsen_mo_exp))';
                                    R_kans_botsen_mo_exp(nanis) = 0;
                                    %                      nanis = find(isnan(R_kans_botsen_pl_exp))';
                                    %retreive background velocity
                                    Vbg_exp = (Mbg.(x_names{x+n_nbg-1}).V);
                                    % no colisions if the velocity is slower than the background velocity
                                    if Vbg_exp(2) > V_local_dist
                                        Nbg_exp(2) = 0;
                                    end
                                    %calculate the amont of collided atoms
                                    R_kans_bots_st(n_nbg) = R_kans_botsen_st_exp(n_nbg)*N_temp_left(n_nbg)*Nbg_exp(1)/sum(Nbg_exp);
                                    R_kans_bots_mo(n_nbg) = R_kans_botsen_mo_exp(n_nbg) *N_temp_left(n_nbg)*Nbg_exp(2)/sum(Nbg_exp);
                                    %cap if no more particles present in the background
                                    if R_kans_bots_mo(n_nbg) > Nbg_exp_N(2)
                                        R_kans_bots_mo(n_nbg) = Nbg_exp_N(2);
                                    end
                                    if R_kans_bots_st(n_nbg) > Nbg_exp_N(1)
                                        R_kans_bots_st(n_nbg) = Nbg_exp_N(1);
                                    end
                                    if Vbg_exp(2) > V_local_dist
                                        R_kans_bots_mo(n_nbg) = 0;
                                    end
                                    
                                    %end of cap if no more particles present in the background
                                    nanis = find(isnan(R_kans_bots_st))';
                                    R_kans_bots_st(nanis) = 0;
                                    nanis = find(isnan(R_kans_bots_mo))';
                                    R_kans_bots_mo(nanis) = 0;
                                    
                                    %calculate the not collided plasma particles
                                    lost = ( R_kans_bots_st(n_nbg)    + R_kans_bots_mo(n_nbg)  ); %+ % sumR_kans_bots_pl:,n_nbg
                                    if lost > N_temp_left(n_nbg)
                                        temp_R_st = R_kans_bots_st(n_nbg)/sum(R_kans_bots_st(n_nbg) + R_kans_bots_mo(n_nbg)) * N_temp_left(n_nbg);
                                        temp_R_mo = R_kans_bots_mo(n_nbg)/sum(R_kans_bots_st(n_nbg) + R_kans_bots_mo(n_nbg)) * N_temp_left(n_nbg);
                                        R_kans_bots_st(n_nbg) = temp_R_st;
                                        R_kans_bots_mo(n_nbg) = temp_R_mo;
                                        lost = N_temp_left(n_nbg);
                                    end
                                    
                                    nanis = find(isnan(lost))';
                                    lost(nanis) = 0;
                                    N_temp_left(n_nbg+1) = N_temp_left(n_nbg)- lost;
                                    
                                    %calculate plasma verlocities after collisions
                                    V_min_mo = ((m_p(P_local+1)-m_p(2))*V_local_dist+2*m_o*(Vbg_exp(1)))/((m_p(P_local+1)+m_p(2)));
                                    V_min_st = ((m_p(P_local+1)-m_p(2))*V_local_dist+2*m_o*Vbg(1))/(m_p(P_local+1)+m_p(2));
                                    V_min_mo = V_local_dist;
                                    
                                    %Round on stepsize
                                    V_min_st = round(V_min_st/(X_delta/T_delta))*(X_delta/T_delta);
                                    V_min_mo = round(V_min_mo/(X_delta/T_delta))*(X_delta/T_delta);
                                   
                                    if Vbg_exp(2) > V_local_dist
                                        V_min_mo = 0;
                                    end
                                    
                                    if V_local_dist > sqrt(5.2*1.6E-19/m_c*2) && V_local_dist < sqrt(E_O_high*2/m_p(1)) && P_local==0
                                        P_local_temp = 2 ;
                                    else
                                        P_local_temp = P_local;
                                    end
                                    
                                    %calculate background verlocities after collisions
                                    Vbg_min_st = ((m_o-m_p(P_local_temp+1))*0+2*m_p(P_local_temp+1)*V_local_dist)/(m_p(P_local_temp+1)+m_o);
                                    Vbg_min_mo(n_nbg) = round(((m_o-m_p(P_local_temp+1))*Vbg_exp(2)+2*m_p(P_local_temp+1)*V_local_dist)/(m_p(P_local_temp+1)+m_o)/(X_delta/T_delta))*(X_delta/T_delta);
                                   
                                    %calculate traveled distance
                                    S_st(n_nbg) = V_local_dist*T_delta*n_nbg/S_exp+V_min_st*T_delta*(S_exp-n_nbg)/S_exp;
                                    S_mo(n_nbg) = V_local_dist*T_delta*n_nbg/S_exp+V_min_mo*T_delta*(S_exp-n_nbg)/S_exp;
                                   
                                    Sbg_st(n_nbg) = V_local_dist*T_delta*n_nbg/S_exp+Vbg_min_st*T_delta*(S_exp-n_nbg)/S_exp;  %aftand waar hij begint + die hij af
                                    Sbg_mo(n_nbg) = Vbg_exp(2)*T_delta*n_nbg/S_exp+Vbg_min_mo(n_nbg)*T_delta*(S_exp-n_nbg)/S_exp;
                                    %round calculated distance to X_delta and set the
                                end %end if x+n_nbg-1 > ((X_end-X_begin)/X_delta)
                                %starting point offset
                                x_new_bots_st = round(S_st/X_delta)*X_delta+X;       % plasma For collisions with the standing Background
                                x_new_bots_mo = round(S_mo/X_delta)*X_delta+X;       % plasma For collisions with the moving Background
                                x_newbg_bots_st = round(Sbg_st/X_delta)*X_delta+X;    % standing background For collisions with the plasma
                                x_newbg_bots_mo = round(Sbg_mo/X_delta)*X_delta+X;    % moving background For collisions with the plasma
                                nanis = find(isnan(x_newbg_bots_mo))';
                                x_newbg_bots_mo(nanis) = 0;
                            end %n_nbg=1:S_exp
                        end %if S_exp==0
                        %amount of particles that dont collide inside the plasma
                        
                        R_kans_nbots = N_local_dist-sum(R_kans_bots_st)-sum(R_kans_bots_mo); %aantal deeltjes dat niet bots
                        S_new_nbots = X+T_delta*V_local_dist; %New position of the not collided
                        %round to stepsize
                        
                        x_new_nbots = round(S_new_nbots/X_delta)*X_delta;   %afronden naar het dichtsbijzijnde X punt. voor niet bots
                        %Write locations into string formation
                        if x_new_nbots > 0
                            X_n_nb = ['X_' myfastint2str(round(x_new_nbots/X_delta+1))];  %for the non collided particles
                        else
                            X_n_nb = ['X_' num2str(round(x_new_nbots/X_delta+1))];  %for the non collided particles
                        end
                          
                        % current position in string format
                        if round(x) == 0
                            Xbg_n_minb = 'X_0';
                        elseif round(x) == 1
                            Xbg_n_minb = 'X_1';
                        else 
                            Xbg_n_minb = ['X_' myfastint2str(round(x))];
                        end
                        
                        %% check moment
                        % %--------------------------------------------------------------------------
                        
                        %% Writing particles into struture M_1 en Mbg
                        %
                        % Non coliding particles
                        %
                        %--------------------------------------------------------------------------
                        %de deeltjes die aankomen
                        if X_TS < cosd(rad-Rad_delta)*x_new_nbots  && X_TS >= cosd(rad-Rad_delta)*X  && R_kans_nbots > 1 && V_local_dist>-1 %&& V_local_dist<8E4
                            Ek_loc = round(V_local_dist/V_minimaal)+1;
                            Ek_arr_N(Ek_loc,P_local+1) = Ek_arr_N(Ek_loc,P_local+1)+R_kans_nbots; %aantal,tijd,type
                        end
                        %de deeltjes die verdwijnen
                        if x_new_nbots >= (X_end-X_delta) || x_new_nbots < 0 || R_kans_nbots < 1 %Check if they disappear
                            
                        else %if they dont disappear
                            %Search if the velocity and number of colisions is already
                            %in the structure
                            
                            i_nb_v = (M_1.(X_n_nb).V == V_local_dist);
                            i_nb_k = (M_1.(X_n_nb).K == K_local);
                            i_nb_p = (M_1.(X_n_nb).P == P_local);
                            if isempty(i_nb_v)==0 %snelehid is niet gevonden in de array
                                loc_nb = find(i_nb_v.*i_nb_k.*i_nb_p == 1) ;%if verlocity and K collisions are equal
                                if loc_nb > 0 %if there is a location
                                    N_loc_nb = M_1.(X_n_nb).N(loc_nb); %Amount of particles on the new location
                                    M_1.(X_n_nb).N(loc_nb) = N_loc_nb+R_kans_nbots; %Add particles to the location
                                else %If there is no existing data make a new entry
                                    %adding the data
                                    i_nb_size = size(M_1.(X_n_nb).V'); %Look at size of the array
                                    M_1.(X_n_nb).V(i_nb_size(1)+1) = V_local_dist; %add velocity
                                    M_1.(X_n_nb).N(i_nb_size(1)+1) = R_kans_nbots; %add amount of particles
                                    M_1.(X_n_nb).K(i_nb_size(1)+1) = K_local; %add value of number of colisions
                                    M_1.(X_n_nb).P(i_nb_size(1)+1) = P_local; %add value of number of colisions
                                end %end loc_nb > 0
                            else %else snelehid is niet gevonden in de array
                                %add if the array has no size
                                i_nb_size = size(M_1.(X_n_nb).V'); %Look at size of the array
                                M_1.(X_n_nb).V(i_nb_size(1)+1) = V_local_dist; %%add velocity
                                M_1.(X_n_nb).N(i_nb_size(1)+1) = R_kans_nbots;  %add amount of particles
                                M_1.(X_n_nb).K(i_nb_size(1)+1) = K_local;%add value of number of colisions
                                M_1.(X_n_nb).P(i_nb_size(1)+1) = P_local;%add value of type of colisions
                            end    %end if empty
                        end  %Check if they disappear
                        for all = 1:S_exp
                            %%
                            %
                            % Colliding particles with a standing background
                            %
                            %--------------------------------------------------------------------------
                            % Xbg_n_minb    %current position
                            % X_n_b_st      % new location plasma
                            % Xbg_n_b_st    % new location backgound
                            %
                            %write into single string and if X positions is 0 change to -0
                            % some values are preset for improving
                            % calucation speed
                            value_X_n_b_st = round(x_new_bots_st(all)/X_delta+1); 
                            if value_X_n_b_st == 0
                                X_n_b_st = 'X_0';
                            elseif value_X_n_b_st == 1
                                X_n_b_st = 'X_1';
                            elseif value_X_n_b_st == ((X_end-X_begin)/X_delta)-2;
                                X_n_b_st = X_str_end2;
                            elseif value_X_n_b_st == ((X_end-X_begin)/X_delta)-1;
                                X_n_b_st = X_str_end1;
                            elseif value_X_n_b_st > 0
                                X_n_b_st = ['X_' myfastint2str(value_X_n_b_st)]; 
                            else
                                X_n_b_st = 'X_-0';
                            end
                            
                            % some values are preset for improving
                            % calucation speed
                            value_X_nbg_b_st_a = round(x_newbg_bots_st(all)/X_delta+1); 
                            if value_X_nbg_b_st_a == 0
                                Xbg_n_b_st = 'X_0';
                            elseif value_X_nbg_b_st_a == 1
                                Xbg_n_b_st = 'X_1';
                            elseif value_X_nbg_b_st_a == ((X_end-X_begin)/X_delta)-2;
                                Xbg_n_b_st = X_str_end2;
                            elseif value_X_nbg_b_st_a == ((X_end-X_begin)/X_delta)-1;
                                Xbg_n_b_st = X_str_end1;
                            elseif value_X_nbg_b_st_a > 0
                                Xbg_n_b_st = ['X_' myfastint2str(value_X_nbg_b_st_a)];
                            else
                                Xbg_n_b_st = 'X_-0';
                            end
                            %de deeltjes die aankomen ( de voorwaarden,
                            %deeltjes die een grotere X_ts hebben dan X_ts
                            % deeltjes die starten vanaf X die kleiners is
                            % dan X_TS
                            % kan geen negatieeve snelheid hebben om te
                            % tellen
                            % meer dan 1 deeltjem gene hele kleine
                            % fracties, scheelt rekentijd.
                            
                            if  X_TS < cosd(rad-Rad_delta)*x_new_bots_st(all) && X_TS >= cosd(rad-Rad_delta)*X && V_local_dist>-1  && R_kans_bots_st(all) > 1
                                Ek_loc = round(V_local_dist/V_minimaal)+1;
                                Ek_arr_N(Ek_loc,P_local+1) = Ek_arr_N(Ek_loc,P_local+1)+R_kans_bots_st(all); %aantal,tijd,type
                            end
                            %de deeltjes die verdwijnen
                            if x_new_bots_st(all) >= (X_end-X_delta)  || x_new_bots_st(all) < 0 || R_kans_bots_st(all) < 1 || X_n_b_st(3)=='-'%Check if they disappear
                                
                            else
                                %Reactions with Oxygen:
                                %---------------------------
                                %NOTE making of CuO
                                %makeing of MO
                                if V_local_dist > sqrt(E_O_low*2/m_p(1)) && V_local_dist < sqrt(E_O_high*2/m_p(1)) && P_local==0
                                    P_local_temp = 2 ;
                                else
                                    P_local_temp = P_local;
                                end
                                %NOTE making of TIO2 MO2
                                if MO2 == 1
                                    if V_local_dist > sqrt(E_O_low*2/m_p(3)) && V_local_dist < sqrt(E_O_high*2/m_p(3)) && P_local==2
                                        P_local_temp = 3 ;
                                    end
                                end
                                %---------------------------
                                %Plocal 0 = M
                                %Plocal 1 = O
                                %P local 2 = MO
                                %P local 3 = MO2
                                %zoek of snelheid al bestaat voor het eventueel bijschijven
                                i_b_v = (M_1.(X_n_b_st).V == V_min_st);
                                i_b_k = (M_1.(X_n_b_st).K == (K_local+1));
                                i_b_p = (M_1.(X_n_b_st).P == P_local_temp);
                                if isempty(i_b_v) == 0 %snelehid is niet gevonden in de array
                                    loc_b=find(i_b_v.*i_b_k.*i_b_p == 1); %vind de locatie waar k en v gelijk zijn
                                    if loc_b > 0 %als er een lcoatie is
                                        N_loc_b = M_1.(X_n_b_st).N(loc_b); %vooraf aantal aanwezige deeltjes bepalen op locatie
                                        M_1.(X_n_b_st).N(loc_b) = N_loc_b+R_kans_bots_st(all); %deeltje bijschrijven op nieuwe locatatie voor niet bots
                                    else
                                        i_b_size = size(M_1.(X_n_b_st).V'); %kijken hoe groot de array is
                                        M_1.(X_n_b_st).V(i_b_size(1)+1) = V_min_st; %nieuwe snelheid bijscrijven
                                        M_1.(X_n_b_st).N(i_b_size(1)+1) = R_kans_bots_st(all); %nieuwe verdeling bijschrijven
                                        M_1.(X_n_b_st).K(i_b_size(1)+1) = K_local+1; %aantal geboste deeltjes
                                        M_1.(X_n_b_st).P(i_b_size(1)+1) = P_local_temp; %type deeltje
                                    end %end loc_b > 0
                                else %geen bestaande plek om data erbij te voegen
                                    %bijschrijven van de nieuwe data
                                    i_b_size = size(M_1.(X_n_b_st).V'); %kijken hoe groot de array is
                                    M_1.(X_n_b_st).V(i_b_size(1)+1) = V_min_st; %nieuwe snelheid bijscrijven
                                    M_1.(X_n_b_st).N(i_b_size(1)+1) = R_kans_bots_st(all); %nieuwe verdeling bijschrijven
                                    M_1.(X_n_b_st).K(i_b_size(1)+1) = K_local+1; %aantal geboste deeltjes
                                    M_1.(X_n_b_st).P(i_b_size(1)+1) = P_local_temp; %aantal geboste deeltjes
                                end %end if empty
                            end %bijschrijven bostingen
                            % clear i_b_v i_b_size i_b_v i_b_k Vbg_new  i_b_p
                            %for the background
                            Xbg_n_minb_curr = ['X_' myfastint2str(round(x+all-1))]; %curret all
                            if x+all-1 > (X_end/X_delta)
                            else
                                if x_newbg_bots_st(all) >= (X_end-X_delta) || x_newbg_bots_st(all) < 0 || Xbg_n_b_st(3)=='-' %Check if they disappear
                                    Mbg.(Xbg_n_minb_curr).N(1) = Mbg.(Xbg_n_minb_curr).N(1)-R_kans_bots_st(all); %remove background particles
                                else %if they move to a new space
                                    %remove background particles
                                    Mbg.(Xbg_n_b_st).N(2) = Mbg.(Xbg_n_b_st).N(2)+R_kans_bots_st(all);
                                    %add particles
                                    Mbg.(Xbg_n_minb_curr).N(1) = Mbg.(Xbg_n_minb_curr).N(1)-R_kans_bots_st(all);
                                    %add new velocities
                                    Vbg_new = (Mbg.(Xbg_n_b_st).V(2)*(Mbg.(Xbg_n_b_st).N(2)-R_kans_bots_st(all))+ R_kans_bots_st(all)*Vbg_min_st)/(Mbg.(Xbg_n_b_st).N(2));
                                    if isnan(Vbg_new) == 1
                                        Vbg_new = 0;
                                    end
                                    Mbg.(Xbg_n_b_st).V(2) = Vbg_new ;
                                end %end vCheck if they disappear
                            end %end if x+all-1 > (X_end-X_delta)
                            %End of the part background
                            %%
                            %
                            % Colliding particles with a moving background
                            %
                            %--------------------------------------------------------------------------
                            % Xbg_n_minb    %huidige positie
                            % X_n_b_mo      % nieuwe locatie van plasma deeltjes die gebotst zijn met stil staand plasma
                            % Xbg_n_b_mo    % nieuwe locatie van background deeltjes die gebotst zijn met stil staand plasma
                            %
                            
                            % Voor de plume
                            % some values are preset for improving
                            % calucation speed
                            value_X_n_b_mo_a = round(x_new_bots_mo(all)/X_delta+1);
                            if value_X_n_b_mo_a == 0
                                X_n_b_mo = 'X_0';
                            elseif value_X_n_b_mo_a == 1
                                X_n_b_mo = 'X_1';
                            elseif value_X_n_b_mo_a == ((X_end-X_begin)/X_delta)-2;
                                X_n_b_mo = X_str_end2;
                            elseif value_X_n_b_mo_a == ((X_end-X_begin)/X_delta)-1;
                                X_n_b_mo = X_str_end1;
                            elseif value_X_n_b_mo_a > 0
                                X_n_b_mo = ['X_' myfastint2str(value_X_n_b_mo_a)];
                            else
                                X_n_b_mo = 'X_-0';
                            end
                            
                            % some values are preset for improving
                            % calucation speed
                            value_X_nbg_b_mo_a = round(x_newbg_bots_mo(all)/X_delta+1);
                            if value_X_nbg_b_mo_a == 0
                                Xbg_n_b_mo = 'X_0';
                            elseif value_X_nbg_b_mo_a == 1   
                                Xbg_n_b_mo = 'X_1';
                            elseif value_X_nbg_b_mo_a == ((X_end-X_begin)/X_delta)-2;
                                Xbg_n_b_mo = X_str_end2;
                            elseif value_X_nbg_b_mo_a == ((X_end-X_begin)/X_delta)-1;
                                Xbg_n_b_mo = X_str_end1;
                            elseif value_X_nbg_b_mo_a > 0
                                Xbg_n_b_mo = ['X_' myfastint2str(value_X_nbg_b_mo_a)];
                            else
                                Xbg_n_b_mo = 'X_-0';
                            end
                                

                            %de deeltjes die aankomen
                            if X_TS < cosd(rad-Rad_delta)*x_new_bots_mo(all)  && X_TS >= cosd(rad-Rad_delta)*X   && R_kans_bots_mo(all) > 1  && V_local_dist>-1 %&& V_local_dist<8E4
                                Ek_loc = round(V_local_dist/V_minimaal)+1;
                                Ek_arr_N(Ek_loc,P_local+1) = Ek_arr_N(Ek_loc,P_local+1)+R_kans_bots_mo(all); %aantal,tijd,type
                            end
                            %de deetljes die verdwijnen
                            if x_new_bots_mo(all) >= (X_end-X_delta) || x_new_bots_mo(all) == X_end || x_new_bots_mo(all) < 0 || R_kans_bots_mo(all) < 1 || X_n_b_mo(3)=='-'%Check if they disappear
                            else
                                %Reactions with Oxygen:
                                %---------------------------
                                %NOTE making of TIO
                                if V_local_dist > sqrt(E_O_low*2/m_p(1)) && V_local_dist < sqrt(E_O_high*2/m_p(1)) && P_local==0
                                    P_local_temp = 2 ;
                                else
                                    P_local_temp = P_local;
                                end
                                %NOTE making of TIO2
                                if MO2 == 1
                                    if V_local_dist > sqrt(E_O_low*2/m_p(3)) && V_local_dist < sqrt(E_O_high*2/m_p(3)) && P_local==2
                                        P_local_temp = 3 ;
                                    end
                                end
                                %---------------------------
                                %Plocal 0 = M
                                %Plocal 1 = O
                                %P local 2 = MO
                                %P local 3 = MO2
                                %same notes as above only in Dutch
                                %zoek of snelheid al bestaat voor het eventueel bijschijven
                                i_b_v = (M_1.(X_n_b_mo).V == V_min_mo);
                                i_b_k = (M_1.(X_n_b_mo).K == (K_local+1));
                                i_b_p = (M_1.(X_n_b_mo).P == P_local_temp);
                                if isempty(i_b_v)==0 %snelehid is niet gevonden in de array
                                    loc_b = find(i_b_v.*i_b_k.*i_b_p == 1); %vind de locatie waar k en v gelijk zijn
                                    if loc_b > 0 %als er een lcoatie is
                                        N_loc_b = M_1.(X_n_b_mo).N(loc_b); %vooraf aantal aanwezige deeltjes bepalen op locatie
                                        M_1.(X_n_b_mo).N(loc_b) = N_loc_b+R_kans_bots_mo(all); %deeltje bijschrijven op nieuwe locatatie voor niet bots
                                    else
                                        i_b_size = size(M_1.(X_n_b_mo).V'); %kijken hoe groot de array is
                                        M_1.(X_n_b_mo).V(i_b_size(1)+1) = V_min_mo; %nieuwe snelheid bijscrijven
                                        M_1.(X_n_b_mo).N(i_b_size(1)+1) = R_kans_bots_mo(all); %nieuwe verdeling bijschrijven
                                        M_1.(X_n_b_mo).K(i_b_size(1)+1) = K_local+1; %aantal geboste deeltjes
                                        M_1.(X_n_b_mo).P(i_b_size(1)+1) = P_local_temp; %type geboste deeltjes
                                    end %end loc_b > 0
                                else %geen bestaande plek om data erbij te voegen
                                    %bijschrijven van de nieuwe data
                                    i_b_size = size(M_1.(X_n_b_mo).V'); %kijken hoe groot de array is
                                    M_1.(X_n_b_mo).V(i_b_size(1)+1) = V_min_mo; %nieuwe snelheid bijscrijven
                                    M_1.(X_n_b_mo).N(i_b_size(1)+1) = R_kans_bots_mo(all); %nieuwe verdeling bijschrijven
                                    M_1.(X_n_b_mo).K(i_b_size(1)+1) = K_local+1; %aantal geboste deeltjes
                                    M_1.(X_n_b_mo).P(i_b_size(1)+1) = P_local_temp; %type geboste deeltjes
                                end %end if empty
                            end %bijschrijven bostingen MO
                            %background
                            if x+all-1 > (X_end/X_delta)
                            else
                                if x_newbg_bots_mo(all) >= (X_end-X_delta) || x_newbg_bots_mo(all) < 0 || Xbg_n_b_mo(3)=='-' %als deeltjes buiten meetgebied komen
                                    Mbg.(Xbg_n_minb_curr).N(2) = Mbg.(Xbg_n_minb_curr).N(2)-R_kans_bots_mo(all); %deeltjes utischrijven
                                else %als deeltjes verplaatsen naar nieuwe ruimte
                                    Mbg.(Xbg_n_b_mo).N(2) = Mbg.(Xbg_n_b_mo).N(2)+R_kans_bots_mo(all);   %bijschrijven achtergrond deeltjes in nieuw vakje vanuit bewegen
                                    Mbg.(Xbg_n_minb_curr).N(2) = Mbg.(Xbg_n_minb_curr).N(2)-R_kans_bots_mo(all); %weghalen deeltjes in huidig vakje
                                    %bijschrijven snelheden
                                    Vbg_new = (Mbg.(Xbg_n_b_mo).V(2)*(Mbg.(Xbg_n_b_mo).N(2)-R_kans_bots_mo(all))+ R_kans_bots_mo(all)*Vbg_min_mo(all))/(Mbg.(Xbg_n_b_mo).N(2));
                                    if isnan(Vbg_new)==1 || Vbg_new == Inf
                                        Vbg_new=0;
                                    end
                                    if Vbg_new > X_end/T_delta
                                        Vbg_new = 600000;
                                    end
                                    Mbg.(Xbg_n_b_mo).V(2) = Vbg_new ;
                                end %end van background als deeltjes uit gebied vliegen
                            end %x+all-1 > (X_end/X_delta)
                            
                        end %end all S_exp
                        
                        
                        % Remove entries smaller then 1
                        if x+all-1 > (X_end/X_delta)
                        else
                            if Mbg.(Xbg_n_minb_curr).N(1) < 1
                                Mbg.(Xbg_n_minb_curr).N(1) = 0;
                            end
                            if Mbg.(Xbg_n_minb_curr).N(2)< 1
                                Mbg.(Xbg_n_minb_curr).N(2) = 0;
                                Mbg.(Xbg_n_minb_curr).V(2) = 0;
                            end
                            if x_newbg_bots_mo(all) >= (X_end-X_delta) || x_newbg_bots_mo(all) < 0 || Xbg_n_b_mo(3)=='-'
                            else
                                if Mbg.(Xbg_n_b_mo).N(2) < 1
                                    Mbg.(Xbg_n_b_mo).N(2) = 0;
                                    Mbg.(Xbg_n_b_mo).V(2) = 0;
                                end
                            end
                        end %x+all-1 > (X_end/X_delta)
                        %  vv;
                    end %end of als V==0 %is Nan
                    %% check moment
                    %--------------------------------------------------------------------------
                    
                end %end of vv
                
                % Move the remaining background particles
                Sbg_nb = Mbg.(Xbg_n_minb).V(2)*T_delta;
                x_newbg_nbots = round(Sbg_nb/X_delta)*X_delta+X;
                if x_newbg_nbots > 0
                    Xbg_n_nb = ['X_' myfastint2str(round(x_newbg_nbots/X_delta+1))];
                else
                    Xbg_n_nb = ['X_' num2str(round(x_newbg_nbots/X_delta+1))];
                end
                if x_newbg_nbots >= (X_end-X_delta) || x_newbg_nbots == X_end || x_newbg_nbots < 0  || Xbg_n_nb(3)=='-' %check if they disappear
                    Mbg.(Xbg_n_minb).N(2) = 0;
                else %is they move to a new space
                    %calculate new velocity
                    Vbg_new = (Mbg.(Xbg_n_nb).N(2)*Mbg.(Xbg_n_nb).V(2) + Mbg.(Xbg_n_minb).V(2) * Mbg.(Xbg_n_minb).N(2))/(Mbg.(Xbg_n_minb).N(2)+Mbg.(Xbg_n_nb).N(2));
                    if isnan(Vbg_new) == 1
                        Vbg_new = 0;
                    end
                    Mbg.(Xbg_n_nb).V(2) = Vbg_new ;
                    Mbg.(Xbg_n_nb).N(2) = Mbg.(Xbg_n_nb).N(2)+Mbg.(Xbg_n_minb).N(2);   %add particles
                    Mbg.(Xbg_n_minb).N(2) = 0; %remove particles
                end
            end %end of fieldnames is small
        end %end of xx
        %% check moment
        
        %
        %% Get data from structure to arrays
        %--------------------------------------------------------------------------
%         clear Nbg_plot_2 Nbg_Vplot_3
%         %pre allocate space
%         N_plot = zeros(round((X_end-X_begin)/X_delta),limbots+1,t_p+t_o-1);
%         V_plot = N_plot;
%         V_plot_1 = N_plot;
%         N_plot_5 = N_plot;
%         X_plot = zeros(round((X_end-X_begin)/X_delta),1)';
%         % 
%         % time, bots, particle
%         %pre allocate space make empty N_plot_5 matrix
%         
%         %get data out of M matrix into array
%         for n = 1:size(x_names)
%             X_plot(n) = (M_1.(x_names{n}).X);
%             Nbg_plot(n) = (Mbg.(x_names{n}).N(1));
%             Nbg_plot_1(n) = (Mbg.(x_names{n}).N(2));
%             Nbg_Vplot_1(n) = (Mbg.(x_names{n}).V(2));
%             for  m = 1:size(M_1.(x_names{n}).V')
%                 k_plot = (M_1.(x_names{n}).K(m));
%                 if k_plot > limbots
%                     k_plot = limbots;
%                 end
%                 p_plot = (M_1.(x_names{n}).P(m));
%                 V_temp = (V_plot(n,k_plot+1,p_plot+1)*N_plot(n,k_plot+1,p_plot+1) + (M_1.(x_names{n}).V(m))*(M_1.(x_names{n}).N(m)))/((M_1.(x_names{n}).N(m))+N_plot(n,k_plot+1,p_plot+1));
%                 N_plot(n,k_plot+1,p_plot+1) = N_plot(n,k_plot+1,p_plot+1)+(M_1.(x_names{n}).N(m));
%                 nanis = find(isnan(V_temp))';
%                 V_temp(nanis) = 0;
%                 V_plot(n,k_plot+1,p_plot+1) = V_temp;
%             end
%         end
%         v_size = size(V_plot) ;
%         N_plot = N_plot(1:v_size(1),1:v_size(2),1:(t_p+t_o-1));
%         
%         % smooth and normlize, 
%         % moving average, standaard function
%         
%         %pre allocate smooth/normalize space
%         N_plot_7 = N_plot.*0;
%         V_plot_7 = N_plot.*0;
%         sizeN5 = size(N_plot);
%         
%         %for all
%         for part=1:sizeN5(3) %for all particles
%             for bots=1:sizeN5(2) %for all bots
%                 if max(N_plot(:,bots,part)) ~= 0 % check if not all values are zero
%                     
%                     % divide into parts
%                     fitplotV = V_plot(:,bots,part);
%                     fitplotN = N_plot(:,bots,part);
%                     fitplotY = fitplotV(fitplotV>0);
%                     fitplotX = X_plot(fitplotV>0);
%                     
%                     % determining smooth value
%                     % for the expectiont of a very low number of distances/spaces occupied
%                     % a lower smooth value is required. 
%                     if sum(fitplotN>1) < 10 %exeption 
%                         smooth_value = 	sum(fitplotN>1)/2;
%                     else
%                         smooth_value = smooth_value_set;
%                     end
%                     
%                     %smooth Velocity
%                     if bots == 1 % for the non collided parts, it is always a linear velocity profile
%                         loc_Velocity = find(fitplotN~=0); %find locations of non-zero elements, only smooth the containers that contain particles
%                         V_fit = polyfit(fitplotX, fitplotY', 1); %a bots=1 always has a linear velocity profile
%                         testV = V_fit(2)+ X_plot(loc_Velocity(1):loc_Velocity(end)) * V_fit(1); % make linear profile
%                         V_plot_7(loc_Velocity(1):loc_Velocity(end),bots,part) = testV; %write to array
%                     else % for non-linear profile, smoothing average
%                         V_plot_7(:,bots,part) = smooth(X_plot,V_plot(:,bots,part),smooth_value+1);
%                     end
%                     
%                     %smooth particles
%                     loc_particles = find(fitplotN~=0); 
%                     %select area of smoothing
%                     if loc_particles(end)+3 < length(fitplotN)
%                         einde_N = loc_particles(end)+3;
%                         Begin_N = loc_particles(1);
%                     else
%                         einde_N = length(fitplotN);
%                         Begin_N = loc_particles(1);
%                     end
%                     
%                     %remove fractions
%                     dd = round(fitplotN(Begin_N:einde_N));
%                     dd(dd < max(dd)/smooth_fraction ) = NaN;
%                     
%                     %fitplotN
%                     testN  = smooth(dd,smooth_value,'moving');
%                     testN(isnan(testN)) = 0;
%                     N_plot_7(Begin_N:einde_N,bots,part) = abs(testN);
%                 end % end of check if not all values are zero
%             end
%         end
%         N_plot_5 = N_plot_7;
%         V_plot_1 = V_plot_7;
%         
%         %----------------------
%         
%         % remove fractions samller then 0
%         N_plot_5(N_plot_5(:,:)<0) = 0;
%         V_plot_1(N_plot_5(:,:)<0) = 0;
%         %rotate matrix
%         Nbg_plot_2 = Nbg_plot_1';
%         Nbg_Vplot_3 = Nbg_Vplot_1';
%         %normalize particles
%         X_plot_5_lim = round((X_TS-X_begin)/X_delta/cosd(rad-Rad_delta)+1);
%         if X_plot_5_lim > (X_end-X_begin)/X_delta
%             X_plot_5_lim = (X_end-X_begin)/X_delta;
%         end
%         for i = 1:v_size(2)
%             for ip = 1:(t_p+t_o-1)
%                 N_plot_5(:,i,ip) = N_plot_5(:,i,ip).*sum(N_plot(1:X_plot_5_lim,i,ip))/sum(N_plot_5(1:X_plot_5_lim,i,ip));
%             end
%         end %end i=1:v_size(2)
%         %normalization factor
%         nor_V = zeros(15,(t_p+t_o-1));
%         V_plot_2 = V_plot_1;
%         for ip = 1:(t_p+t_o-1)
%             nor_V_temp = sum(N_plot(:,:,ip).*V_plot(:,:,ip))./sum(N_plot_5(:,:,ip).*V_plot_2(:,:,ip));
%             size_norv = size(nor_V_temp);
%             nor_V(1:size_norv(2),ip) = nor_V_temp;
%         end
%         nanis = find(isnan(nor_V))';
%         nor_V(nanis) = 0;
%         %normalize velocity
%         for i = 1:v_size(2)
%             for ip = 1:(t_p+t_o-1)
%                 V_plot_1(:,i,ip) = V_plot_2(:,i,ip)*nor_V(i,ip);
%             end
%         end
%         %solve dividing to zero and 1 over 0
%         nanis = find(isnan(Nbg_plot_2))';
%         Nbg_plot_2(nanis) = 0;
%         nanis = find(isnan(N_plot_5))';
%         N_plot_5(nanis) = 0;
%         nanis = find(isnan(V_plot_1))';
%         V_plot_1(nanis) = 0;
%         nanis = find(isnan(Nbg_Vplot_3))';
%         Nbg_Vplot_3(nanis) = 0;
%         %sum data points
%         Nbg_size = size(Nbg_plot');
%         Nbg_plot_som = ((Nbg_plot_2(1:Nbg_size(1),1)+Nbg_plot'));
%         %remove smallen than 1E1
%         N_plot_5(N_plot_5(:,:,:)<1E1) = 0;
%         V_plot_1(N_plot_5(:,:,:)<1E1) = 0;
%         
%         clear M_1
%         %% Calculate physical properties
%         %--------------------------------------------------------------------------
%         %count  Plasma atoms
%         N_diff_p(t) = N_atoom-sum(sum(sum(N_plot_5)));
%         Np(t) = sum(sum(sum(N_plot_5)));
%         %count difference of background atoms
%         N_diff_bg(t) = Total_O2-sum(Nbg_plot_som);
%         Nbg_st(t) = sum(sum(Nbg_plot));
%         Nbg_mo(t) = (sum(Nbg_plot_2(:,1)));
%         %calculate total kinectic energy for the plasma
%         Ek_plasma(t) = sum(sum(sum(N_plot_5.*V_plot_1.^2)))*0.5*m_c;
%         Ek_bg(t) = sum(sum(Nbg_plot_2.*Nbg_Vplot_3.^2))*0.5*m_c;
%         %calculate impuls
%         p_plasma(t) = sum(sum(sum(N_plot_5.*V_plot_1))).*m_c;
%         p_bg(t) = sum(Nbg_plot_2.*Nbg_Vplot_3).*m_o;
%         %% Write new data into M
%         %--------------------------------------------------------------------------
%         for n = 1:size(x_names) %for all X
%             M_1.(x_names{n}).X = M.(x_names{n}).X;
%             for o = 1:v_size(2) %for all k
%                 for oo = 1:(t_p+t_o-1)
%                     if  N_plot_5( n,o,oo)<1
%                     else
%                         ssx_2 = size(fieldnames(M_1.(x_names{n})));
%                         if  ssx_2(1) <= 1
%                             s_size(1)=1;
%                         else
%                             s_size=size(M_1.(x_names{n}).N');
%                         end
%                         M_1.(x_names{n}).N(s_size(1)+1) = N_plot_5( n,o,oo);
%                         M_1.(x_names{n}).V(s_size(1)+1) = V_plot_1(n,o,oo);
%                         M_1.(x_names{n}).K(s_size(1)+1) = o-1;
%                         M_1.(x_names{n}).P(s_size(1)+1) = oo-1;
%                     end
%                 end
%             end
%             Mbg.(x_names{n}).N(1) = Nbg_plot(n);
%             Mbg.(x_names{n}).N(2) = Nbg_plot_2(n);
%             Mbg.(x_names{n}).V = [0 Nbg_Vplot_3(n)]; %was eerst 2
%             Mbg.(x_names{n}).K = 0;
%             Mbg.(x_names{n}).VOL = Mbg.(x_names{n}).VOL;
%         end
%         % clear data for new loop cycle
%         clear M
%         M = M_1;
%         clear M_1
%         %size of the matrix
%         Max_x_size = size(X_plot');
%         %total plasma
%         Total = 0;
%         for ttg = 1:(t_p+t_o-1)
%             Total = Total+sum(N_plot_5(:,:,ttg)');
%         end
%         stap(t) = t;
%         t               %Display the finished timestep
%         names = fieldnames(M);
%         if sum(Total) > 0
%             k_numb = size(N_plot_5);
%             %determine maxima
%             [N_max_t(t) N_loc_t(t)] = max(smooth(Total));
%             X_max_t(t)=X_plot(N_loc_t(t));
%             [N_max_0(t) N_loc_0(t)] = max(smooth(N_plot_5(:,1)));
%             X_max_0(t) = X_plot(N_loc_0(t));
%             V_X_max_0(t) = V_plot_1(N_loc_0(t),1);
%             if k_numb(2) >= 2
%                 [N_max_1(t) N_loc_1(t)] = max(smooth(N_plot_5(:,2)));
%                 X_max_1(t) = X_plot(N_loc_1(t));
%                 V_X_max_1(t) = V_plot_1(N_loc_1(t),2);
%             else
%                 X_max_1(t) = 0 ;
%                 V_X_max_1(t) = 0 ;
%             end
%             if k_numb(2) >= 3
%                 [N_max_2(t) N_loc_2(t)]=max(smooth(N_plot_5(:,3)));
%                 X_max_2(t)=X_plot(N_loc_2(t));
%                 V_X_max_2(t)=V_plot_1(N_loc_2(t),3);
%             else
%                 X_max_2(t) = 0;
%                 V_X_max_2(t) = 0 ;
%             end
%             if k_numb(2) >= 4
%                 [N_max_2(t) N_loc_3(t)]=max(smooth(N_plot_5(:,4)));
%                 X_max_3(t)=X_plot(N_loc_3(t));
%                 V_X_max_3(t)=V_plot_1(N_loc_3(t),4);
%             else
%                 X_max_3(t) = 0 ;
%                 V_X_max_3(t) = 0 ;
%             end
%             if k_numb(2) >= 5
%                 [N_max_2(t) N_loc_4(t)]=max(smooth(N_plot_5(:,5)));
%                 X_max_4(t) = X_plot(N_loc_4(t));
%                 V_X_max_4(t) = V_plot_1(N_loc_4(t),5);
%             else
%                 X_max_4(t) = 0 ;
%                 V_X_max_4(t) = 0 ;
%             end
%             %determine front
%             %all parnum2strles
%             Total_front(t) = max((Total~=0).*(1:length(Total))); %front pixel
%             Total_front_1(t) = max(((Total)>max(Total)/10).*(1:length(Total)));  %10%max
%             Total_front_2(t) = max(((Total)>max(Total)/5).*(1:length(Total)));  %10%max
%             Total_front_3(t) = max(((Total)>max(Total)/3).*(1:length(Total)));  %10%max
%             Total_front_4(t) = max(((Total)>max(Total)/2.5).*(1:length(Total)));  %10%max
%             %Metal particles
%             Total_M = sum(N_plot_5(:,:,1),2)';
%             Total_front_M(t) = max((Total_M~=0).*(1:length(Total_M))); %front
%             Total_front_M1(t) = max(((Total_M)>max(Total_M)/10).*(1:length(Total_M)));  %10%max
%             Total_front_M2(t) = max(((Total_M)>max(Total_M)/5).*(1:length(Total_M)));  %10%max
%             Total_front_M3(t) = max(((Total_M)>max(Total_M)/3).*(1:length(Total_M)));  %10%max
%             Total_front_M4(t) = max(((Total_M)>max(Total_M)/2.5).*(1:length(Total_M)));  %10%max
%             %metal particle
%             [ym yl] = max(Total_M);
%             if max(find(Total_M > ym/10 )) > 0
%                 Total_front_M_2(t) = max(find(Total_M > ym/10 )) ;
%             else
%                 Total_front_M_2(t) = 0;
%             end
%             %Oxide particles
%             Total_O = sum(N_plot_5(:,:,2),2)';
%             Total_front_O(t) = max((Total_O~=0).*(1:length(Total_O))); %front
%             Total_front_O1(t) = max(((Total_O)>max(Total_O)/10).*(1:length(Total_O)));  %10%max
%             Total_front_O2(t) = max(((Total_O)>max(Total_O)/5).*(1:length(Total_O)));  %10%max
%             Total_front_O3(t) = max(((Total_O)>max(Total_O)/3).*(1:length(Total_O)));  %10%max
%             Total_front_O4(t) = max(((Total_O)>max(Total_O)/2.5).*(1:length(Total_O)));  %10%max
%             [ym yl]= max(Total_O);
%             if max(find(Total_O > ym/10 )) > 0
%                 Total_front_O_2(t) = max(find(Total_O > ym/10 )) ;
%             else
%                 Total_front_O_2(t) = 0;
%             end
%             
%             for nn = 1:size(Total')
%                 temp(nn) = sum(Total(1:nn));
%                 if temp(nn) > N_atoom*0.8;
%                     Total_front_deel_2(t) = nn;
%                     break
%                 else
%                     Total_front_deel_2(t) = X_end/X_delta;
%                 end
%             end
%             for nn = 1:size(Total')
%                 temp(nn) = sum(Total(1:nn));
%                 if temp(nn) > N_atoom*0.7;
%                     Total_front_deel_3(t) = nn;
%                     break
%                 else
%                     Total_front_deel_3(t) = X_end/X_delta;
%                 end
%             end
%             for nn = 1:size(Total')
%                 temp(nn) = sum(Total(1:nn));
%                 if temp(nn) > N_atoom*0.9;
%                     Total_front_deel_1(t) = nn;
%                     break
%                 else
%                     Total_front_deel_1(t) = X_end/X_delta;
%                 end
%             end
%             for nn = 1:size(Total')
%                 temp(nn) = sum(Total(1:nn));
%                 if  temp(nn) > N_atoom*0.6;
%                     Total_front_deel_4(t) = nn;
%                     break
%                 else
%                     Total_front_deel_4(t) = X_end/X_delta;
%                 end
%             end
%             %make 1 array for all t
%             Totaal_t(:,t)=Total;
%             if t>1
%                 if sum(Totaal_t(:,t-1)) < sum(Totaal_t(:,t)) ;
%                 end
%                 sum(Totaal_t(:,t-1)) - sum(Totaal_t(:,t))  ;
%             end
%             %% make plots
%             %only plot the center propagation line
%             if rad == Rad_delta && make_1d_plots == 1
%                 % Plot N/V versus distance
%                 %--------------------------------------------------------------------------
%                 fig = figure;
%                 semilogy(X_plot,V_plot_1(:,:,1),'kx--','LineWidth',1)
%                 hold on
%                 semilogy(X_plot,Total(1:Max_x_size)./Vol_T,'cx-','LineWidth',1)
%                 semilogy(X_plot,N_plot_5(1:Max_x_size,1)./Vol_T','bx-','LineWidth',1)
%                 semilogy(X_plot,N_plot_5(1:Max_x_size ,2)./Vol_T','rx-','LineWidth',1)
%                 if k_numb(2) >= 3
%                     semilogy(X_plot,N_plot_5(1:Max_x_size ,3)./Vol_T','gx-','LineWidth',1)
%                 end
%                 if k_numb(2) >= 4
%                     semilogy(X_plot,N_plot_5(1:Max_x_size ,4)./Vol_T','mx-','LineWidth',1)
%                 end
%                 if k_numb(2) >= 5
%                     semilogy(X_plot,N_plot_5(1:Max_x_size ,5)./Vol_T','yx-','LineWidth',1)
%                 end
%                 hold off
%                 xlabel('Distance [m]')
%                 ylabel('Density')
%                 xlim([0 0.05])
%                 title(['time ' num2str(t*T_delta) ' s'])
%                 ylim([1E15 1E25])
%                 if t < 10
%                     legend( 'plasma' ,'0 bots', '1 bots' ,'2 bots','3 bots','4 bots','Location','southeast')
%                 end  %end if t < 10
%                 set(gca,'FontSize',16)
%                 set(gca,'LineWidth',2)
%                 set(get(gca,'xlabel'),'FontSize',16)
%                 set(get(gca,'ylabel'),'FontSize',16)
%                 set(get(gca,'title'),'FontSize',16)
%                 set(findobj('Type','line'),'MarkerSize',10)
%                 set(findobj('Type','line'),'LineWidth',2)
%                 if rad < 10
%                     legend('plasma' ,'0 bots', '1 bots' ,'2 bots','3 bots','4 bots','Location','southeast')
%                     saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\density_log\T0' tijd_string '.jpg'])
%                     saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\density_log\T0' tijd_string '.fig'])
%                 else
%                     saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\density_log\T' tijd_string '.jpg'])
%                     saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\density_log\T' tijd_string '.fig'])
%                 end
%                 close
%                 % Plot N versus distance
%                 %--------------------------------------------------------------------------
%                 fig = figure;
%                 plot(X_plot(Total>1),Total(Total>1),'cx-','LineWidth',1)
%                 if show_BG_inplot == 1
%                     plot(X_plot,Nbg_plot_som,'kx--','LineWidth',1)
%                 end
%                 hold on
%                 plot(X_plot(N_plot_5(:,1)>1),N_plot_5(N_plot_5(:,1)>1,1),'bx-','LineWidth',1)
%                 plot(X_plot(N_plot_5(:,2)>1),N_plot_5(N_plot_5(:,2)>1,2),'rx-','LineWidth',1)
%                 if k_numb(2) >= 3
%                     plot(X_plot(N_plot_5(:,3)>1),N_plot_5(N_plot_5(:,3)>1,3),'gx-','LineWidth',1)
%                 end
%                 if k_numb(2) >= 4
%                     plot(X_plot(N_plot_5(:,4)>1),N_plot_5(N_plot_5(:,4)>1,4),'mx-','LineWidth',1)
%                 end
%                 if k_numb(2) >= 5
%                     plot(X_plot(N_plot_5(:,5)>1),N_plot_5(N_plot_5(:,5)>1 ,5),'yx-','LineWidth',1)
%                 end
%                 hold off
%                 xlabel('Distance [m]')
%                 ylabel('Particles')
%                 xlim([0 0.05])
%                 title(['time ' num2str(t*T_delta) ' s'])
%                 ylim([0 5E12])
%                 if t < 10
%                     legend('plasma' ,'0 bots', '1 bots' ,'2 bots','3 bots','4 bots','Location','southeast')
%                 end  %end if t < 10
%                 set(gca,'FontSize',16)
%                 set(gca,'LineWidth',2)
%                 set(get(gca,'xlabel'),'FontSize',16)
%                 set(get(gca,'ylabel'),'FontSize',16)
%                 set(get(gca,'title'),'FontSize',16)
%                 set(findobj('Type','line'),'MarkerSize',10)
%                 set(findobj('Type','line'),'LineWidth',2)
%                 if rad < 10
%                     saveas(fig,[save_loc folder '\t_p_0\Rad0'  rad_string '\img\particles\T0' tijd_string '.jpg'])
%                     saveas(fig,[save_loc folder '\t_p_0\Rad0'  rad_string '\img\particles\T0' tijd_string '.fig'])
%                 else
%                     saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\particles\T' tijd_string '.jpg'])
%                     saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\particles\T' tijd_string '.fig'])
%                 end
%                 close
%                 
%                 % Plot N versus distance LOG %sum all
%                 %--------------------------------------------------------------------------
%                 fig = figure;
%                 semilogy(X_plot,Total,'cx-','LineWidth',1)
%                 hold on
%                 if show_BG_inplot == 1
%                     semilogy(X_plot,Nbg_plot_som,'kx--','LineWidth',1)
%                 end
%                 semilogy(X_plot,sum(N_plot_5(:,1,:),3),'bx-','LineWidth',1)
%                 semilogy(X_plot,sum(N_plot_5(:,2,:),3),'rx-','LineWidth',1)
%                 if k_numb(2) >= 3
%                     semilogy(X_plot,sum(N_plot_5(:,3,:),3),'gx-','LineWidth',1)
%                 end
%                 if k_numb(2) >= 4
%                     semilogy(X_plot,sum(N_plot_5(:,4,:),3),'mx-','LineWidth',1)
%                 end
%                 if k_numb(2) >= 5
%                     semilogy(X_plot,sum(N_plot_5(:,5,:),3),'yx-','LineWidth',1)
%                 end
%                 hold off
%                 xlabel('Distance X [m]')
%                 ylabel('Particles')
%                 xlim([0 0.05])
%                 title(['time ' num2str(t*T_delta) ' s Particle All'])
%                 ylim([1E7 5E13])
%                 if t < 10
%                     legend('plasma' ,'0 bots', '1 bots' ,'2 bots','3 bots','4 bots','Location','southeast')
%                 end  %end if t < 10
%                 set(gca,'FontSize',16)
%                 set(gca,'LineWidth',2)
%                 set(get(gca,'xlabel'),'FontSize',16)
%                 set(get(gca,'ylabel'),'FontSize',16)
%                 set(get(gca,'title'),'FontSize',16)
%                 set(findobj('Type','line'),'MarkerSize',10)
%                 set(findobj('Type','line'),'LineWidth',2)
%                 if rad < 10
%                     saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\particles_log\T0' tijd_string '.jpg'])
%                     saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\particles_log\T0' tijd_string '.fig'])
%                 else
%                     saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\particles_log\T' tijd_string '.jpg'])
%                     saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\particles_log\T' tijd_string '.fig'])
%                 end
%                 close
%                 
%                  % Plot N versus distance LOG %sum all background density
%                 %--------------------------------------------------------------------------
%                 fig  =figure;
%                 for x = 1:round((X_end-X_begin)/X_delta);
%                     Vol_xy(x) = 4 * pi/3 * (((x)*X_delta)^3-((x-1)*X_delta)^3) * (-cosd(rad)+cosd(rad-Rad_delta));
%                 end
%                 
%                 plot(X_plot,Total./Vol_xy,'cx-','LineWidth',1)
%                 hold on
%                 if show_BG_inplot == 1
%                     plot(X_plot,Nbg_plot_som'./Vol_xy,'kx--','LineWidth',1)
%                 end
%                 plot(X_plot,sum(N_plot_5(:,1,:),3)'./Vol_xy,'bx-','LineWidth',1)
%                 plot(X_plot,sum(N_plot_5(:,2,:),3)'./Vol_xy,'rx-','LineWidth',1)
%                 if k_numb(2) >= 3
%                     plot(X_plot,sum(N_plot_5(:,3,:),3)'./Vol_xy,'gx-','LineWidth',1)
%                 end
%                 if k_numb(2) >= 4
%                     plot(X_plot,sum(N_plot_5(:,4,:),3)'./Vol_xy,'mx-','LineWidth',1)
%                 end
%                 if k_numb(2) >= 5
%                     plot(X_plot,sum(N_plot_5(:,5,:),3)'./Vol_xy,'yx-','LineWidth',1)
%                 end
%                 hold off
%                 xlabel('Distance X [m]')
%                 ylabel('density')
%                 xlim([0 0.05])
%                 title(['time ' num2str(t*T_delta) ' density'])
%                 if t < 10
%                     legend('plasma' ,'0 bots', '1 bots' ,'2 bots','3 bots','4 bots','Location','southeast')
%                 end  %end if t < 10
%                 set(gca,'FontSize',16)
%                 set(gca,'LineWidth',2)
%                 set(get(gca,'xlabel'),'FontSize',16)
%                 set(get(gca,'ylabel'),'FontSize',16)
%                 set(get(gca,'title'),'FontSize',16)
%                 set(findobj('Type','line'),'MarkerSize',10)
%                 set(findobj('Type','line'),'LineWidth',2)
%                 if rad < 10
%                     saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\density\T0' tijd_string '.jpg'])
%                     saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\density\T0' tijd_string '.fig'])
%                 else
%                     saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\density\T' tijd_string '.jpg'])
%                     saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\density\T' tijd_string '.fig'])
%                 end
%                 close
%                 
%                 
%                 
%                 % Plot N versus distance O2
%                 %--------------------------------------------------------------------------
%                 if t_p>1
%                     fig = figure;
%                     semilogy(X_plot,sum(sum(N_plot_5(:,:,2),3),2),'cx-','LineWidth',1)
%                     hold on
%                     if show_BG_inplot == 1
%                         semilogy(X_plot,Nbg_plot_som,'kx--','LineWidth',1)
%                     end
%                     semilogy(X_plot,N_plot_5(:,1,2),'bx-','LineWidth',1)
%                     semilogy(X_plot,N_plot_5(:,2,2),'rx-','LineWidth',1)
%                     if k_numb(2) >= 3
%                         semilogy(X_plot,N_plot_5(:,3,2),'gx-','LineWidth',1)
%                     end
%                     if k_numb(2) >= 4
%                         semilogy(X_plot,N_plot_5(:,4,2),'mx-','LineWidth',1)
%                     end
%                     if k_numb(2) >= 5
%                         semilogy(X_plot,N_plot_5(:,5,2),'yx-','LineWidth',1)
%                     end
%                     hold off
%                     xlabel('Distance [m]')
%                     ylabel('Particles')
%                     xlim([0 0.05])
%                     title(['time ' num2str(t*T_delta) ' s Particle 2'])
%                     ylim([1E7 5E13])
%                     set(gca,'FontSize',16)
%                     set(gca,'LineWidth',2)
%                     set(get(gca,'xlabel'),'FontSize',16)
%                     set(get(gca,'ylabel'),'FontSize',16)
%                     set(get(gca,'title'),'FontSize',16)
%                     set(findobj('Type','line'),'MarkerSize',10)
%                     set(findobj('Type','line'),'LineWidth',2)
%                     if t < 10
%                         legend('plasma' ,'0 bots', '1 bots' ,'2 bots','3 bots','4 bots','Location','southeast')
%                     end  %end if t < 10
%                     if rad < 10
%                         saveas(fig,[save_loc folder '\t_p_' myfastint2str(2) '\Rad0' rad_string '\img\particles_log_O\T0' tijd_string '.jpg'])
%                         saveas(fig,[save_loc folder '\t_p_' myfastint2str(2) '\Rad0' rad_string '\img\particles_log_O\T0' tijd_string '.fig'])
%                     else
%                         saveas(fig,[save_loc folder '\t_p_' myfastint2str(2) '\Rad' rad_string '\img\particles_log_O\T' tijd_string '.jpg'])
%                         saveas(fig,[save_loc folder '\t_p_' myfastint2str(2) '\Rad' rad_string '\img\particles_log_O\T' tijd_string '.fig'])
%                     end
%                     close
%                     
%                 end%t_p>1
%                 % Plot N versus distance Cu
%                 %--------------------------------------------------------------------------
%                 if t_p>1
%                     fig = figure;
%                     semilogy(X_plot,sum(sum(N_plot_5(:,:,1),3),2),'cx-','LineWidth',1)
%                     hold on
%                     if show_BG_inplot == 1
%                         semilogy(X_plot,Nbg_plot_som,'kx--','LineWidth',1)
%                     end
%                     semilogy(X_plot,N_plot_5(:,1,1),'bx-','LineWidth',1)
%                     semilogy(X_plot,N_plot_5(:,2,1),'rx-','LineWidth',1)
%                     if k_numb(2) >= 3
%                         semilogy(X_plot,N_plot_5(:,3,1),'gx-','LineWidth',1)
%                     end
%                     if k_numb(2) >= 4
%                         semilogy(X_plot,N_plot_5(:,4,1),'mx-','LineWidth',1)
%                     end
%                     if k_numb(2) >= 5
%                         semilogy(X_plot,N_plot_5(:,5,1),'yx-','LineWidth',1)
%                     end
%                     hold off
%                     xlabel('Distance [m]')
%                     ylabel('Particles')
%                     xlim([0 0.05])
%                     title(['time ' num2str(t*T_delta) ' s Particle 1'])
%                     ylim([1E7 5E13])
%                     set(gca,'FontSize',16)
%                     set(gca,'LineWidth',2)
%                     set(get(gca,'xlabel'),'FontSize',16)
%                     set(get(gca,'ylabel'),'FontSize',16)
%                     set(get(gca,'title'),'FontSize',16)
%                     set(findobj('Type','line'),'MarkerSize',10)
%                     set(findobj('Type','line'),'LineWidth',2)
%                     if t < 10
%                         legend('plasma' ,'0 bots', '1 bots' ,'2 bots','3 bots','4 bots','Location','southeast')
%                     end  %end if t < 10
%                     if rad < 10
%                         saveas(fig,[save_loc folder '\t_p_' myfastint2str(1) '\Rad0' rad_string '\img\particles_log_P\T0' tijd_string '.jpg'])
%                         saveas(fig,[save_loc folder '\t_p_' myfastint2str(1) '\Rad0' rad_string '\img\particles_log_P\T0' tijd_string '.fig'])
%                     else
%                         saveas(fig,[save_loc folder '\t_p_' myfastint2str(1) '\Rad' rad_string '\img\particles_log_P\T' tijd_string '.jpg'])
%                         saveas(fig,[save_loc folder '\t_p_' myfastint2str(1) '\Rad' rad_string '\img\particles_log_P\T' tijd_string '.fig'])
%                     end
%                     close
%                     
%                 end%t_p>1
%                 
%                 % Plot N versus distance Cu nonlog
%                 %--------------------------------------------------------------------------
%                 if t_p>1
%                     fig=figure;
%                     
%                     plot(X_plot,sum(sum(N_plot_5(:,:,1),3),2),'cx-','LineWidth',1)
%                     hold on
%                     if show_BG_inplot == 1
%                         semilogy(X_plot,Nbg_plot_som,'kx--','LineWidth',1)
%                     end
%                     plot(X_plot,N_plot_5(:,1,1),'bx-','LineWidth',1)
%                     plot(X_plot,N_plot_5(:,2,1),'rx-','LineWidth',1)
%                     if k_numb(2) >= 3
%                         plot(X_plot,N_plot_5(:,3,1),'gx-','LineWidth',1)
%                     end
%                     if k_numb(2) >= 4
%                         plot(X_plot,N_plot_5(:,4,1),'mx-','LineWidth',1)
%                     end
%                     if k_numb(2) >= 5
%                         plot(X_plot,N_plot_5(:,5,1),'yx-','LineWidth',1)
%                     end
%                     hold off
%                     xlabel('Distance [m]')
%                     ylabel('Particles')
%                     xlim([0 0.05])
%                     title(['time ' num2str(t*T_delta) ' s Particle 1'])
%                     ylim([1E7 1E13])
%                     set(gca,'FontSize',16)
%                     set(gca,'LineWidth',2)
%                     set(get(gca,'xlabel'),'FontSize',16)
%                     set(get(gca,'ylabel'),'FontSize',16)
%                     set(get(gca,'title'),'FontSize',16)
%                     set(findobj('Type','line'),'MarkerSize',10)
%                     set(findobj('Type','line'),'LineWidth',2)
%                     if t < 10
%                         legend('plasma' ,'0 bots', '1 bots' ,'2 bots','3 bots','4 bots','Location','southeast')
%                     end  %end if t < 10
%                     if rad < 10
%                         saveas(fig,[save_loc folder '\t_p_' myfastint2str(1) '\Rad0' rad_string '\img\particles\T0' tijd_string '.jpg'])
%                         saveas(fig,[save_loc folder '\t_p_' myfastint2str(1) '\Rad0' rad_string '\img\particles\T0' tijd_string '.fig'])
%                     else
%                         saveas(fig,[save_loc folder '\t_p_' myfastint2str(1) '\Rad' rad_string '\img\particles\T' tijd_string '.jpg'])
%                         saveas(fig,[save_loc folder '\t_p_' myfastint2str(1) '\Rad' rad_string '\img\particles\T' tijd_string '.fig'])
%                     end
%                     close
%                 end%t_p>1
%                 
%                 %Plot N versus distance CuO
%                 %--------------------------------------------------------------------------
%                 if t_p>2
%                     fig=figure;
%                     
%                     semilogy(X_plot,sum(sum(N_plot_5(:,:,3),3),2),'cx-','LineWidth',1)
%                     hold on
%                     if show_BG_inplot == 1
%                         semilogy(X_plot,Nbg_plot_som,'kx--','LineWidth',1)
%                     end
%                     semilogy(X_plot,N_plot_5(:,1,3),'bx-','LineWidth',1)
%                     semilogy(X_plot,N_plot_5(:,2,3),'rx-','LineWidth',1)
%                     if k_numb(2) >= 3
%                         semilogy(X_plot,N_plot_5(:,3,3),'gx-','LineWidth',1)
%                     end
%                     if k_numb(2) >= 4
%                         semilogy(X_plot,N_plot_5(:,4,3),'mx-','LineWidth',1)
%                     end
%                     if k_numb(2) >= 5
%                         semilogy(X_plot,N_plot_5(:,5,3),'yx-','LineWidth',1)
%                     end
%                     hold off
%                     xlabel('Distance [m]')
%                     ylabel('Particles')
%                     xlim([0 0.05])
%                     title(['time ' num2str(t*T_delta) ' s Particle 3'])
%                     ylim([1E7 5E13])
%                     set(gca,'FontSize',16)
%                     set(gca,'LineWidth',2)
%                     set(get(gca,'xlabel'),'FontSize',16)
%                     set(get(gca,'ylabel'),'FontSize',16)
%                     set(get(gca,'title'),'FontSize',16)
%                     set(findobj('Type','line'),'MarkerSize',10)
%                     set(findobj('Type','line'),'LineWidth',2)
%                     if t < 10
%                         legend('plasma' ,'0 bots', '1 bots' ,'2 bots','3 bots','4 bots','Location','southeast')
%                     end  %end if t < 10
%                     if rad < 10
%                         saveas(fig,[save_loc folder '\t_p_' myfastint2str(3) '\Rad0' rad_string '\img\particles_log_CuO\T0' tijd_string '.jpg'])
%                         saveas(fig,[save_loc folder '\t_p_' myfastint2str(3) '\Rad0' rad_string '\img\particles_log_CuO\T0' tijd_string '.fig'])
%                     else
%                         saveas(fig,[save_loc folder '\t_p_' myfastint2str(3) '\Rad' rad_string '\img\particles_log_CuO\T' tijd_string '.jpg'])
%                         saveas(fig,[save_loc folder '\t_p_' myfastint2str(3) '\Rad' rad_string '\img\particles_log_CuO\T' tijd_string '.fig'])
%                     end
%                     close
%                 end%t_p>2
%                 
%                 % Plot N versus distance MO2
%                 %--------------------------------------------------------------------------
%                 if t_p>=3
%                     fig=figure;
%                     
%                     semilogy(X_plot,sum(sum(N_plot_5(:,:,4),4),2),'cx-','LineWidth',1)
%                     hold on
%                     if show_BG_inplot == 1 
%                         semilogy(X_plot,Nbg_plot_som,'kx--','LineWidth',1)
%                     end
%                     semilogy(X_plot,N_plot_5(:,1,4),'bx-','LineWidth',1)
%                     semilogy(X_plot,N_plot_5(:,2,4),'rx-','LineWidth',1)
%                     if k_numb(2) >= 3
%                         semilogy(X_plot,N_plot_5(:,3,4),'gx-','LineWidth',1)
%                     end
%                     if k_numb(2) >= 4
%                         semilogy(X_plot,N_plot_5(:,4,4),'mx-','LineWidth',1)
%                     end
%                     if k_numb(2) >= 5
%                         semilogy(X_plot,N_plot_5(:,5,4),'yx-','LineWidth',1)
%                     end
%                     hold off
%                     xlabel('Distance [m]')
%                     ylabel('Particles')
%                     xlim([0 0.05])
%                     title(['time ' num2str(t*T_delta) ' s Particle 4'])
%                     ylim([1E7 5E13])
%                     set(gca,'FontSize',16)
%                     set(gca,'LineWidth',2)
%                     set(get(gca,'xlabel'),'FontSize',16)
%                     set(get(gca,'ylabel'),'FontSize',16)
%                     set(get(gca,'title'),'FontSize',16)
%                     set(findobj('Type','line'),'MarkerSize',10)
%                     set(findobj('Type','line'),'LineWidth',2)
%                     if t < 10
%                         legend('plasma' ,'0 bots', '1 bots' ,'2 bots','3 bots','4 bots','Location','southeast')
%                     end  %end if t < 10
%                     if rad < 10
%                         saveas(fig,[save_loc folder '\t_p_' myfastint2str(4) '\Rad0' rad_string '\img\particles_log_CuO\T0' tijd_string '.jpg'])
%                         saveas(fig,[save_loc folder '\t_p_' myfastint2str(4) '\Rad0' rad_string '\img\particles_log_CuO\T0' tijd_string '.fig'])
%                     else
%                         saveas(fig,[save_loc folder '\t_p_' myfastint2str(4) '\Rad' rad_string '\img\particles_log_CuO\T' tijd_string '.jpg'])
%                         saveas(fig,[save_loc folder '\t_p_' myfastint2str(4) '\Rad' rad_string '\img\particles_log_CuO\T' tijd_string '.fig'])
%                     end
%                     close
%                 end%t_p>3
%             end %if (rad/Rad_delta)==1
%         end %end if total is kleiner dan 0
%         %end %end if only afbeeldingen op bepaalde tijdstappen
%         %write individual values for collisions
%         for gh=1:(t_p+t_o-1)
%             if rad < 10
%                 mkdir([save_loc folder '\t_p_' myfastint2str(gh)  '\Rad0' rad_string '\N'])
%                 mkdir([save_loc folder '\t_p_' myfastint2str(gh)  '\Rad0' rad_string '\V'])
%             else
%                 mkdir([save_loc folder '\t_p_' myfastint2str(gh)  '\Rad' rad_string '\N'])
%                 mkdir([save_loc folder '\t_p_' myfastint2str(gh)  '\Rad' rad_string '\V'])
%             end
%             %wrting down the particles
%             if rad < 10
%                 fid=fopen([save_loc folder '\t_p_' myfastint2str(gh)  '\Rad0' rad_string '\N\tijd_' tijd_string '.txt'],'w');
%                 fclose(fid);
%                 csvwrite([save_loc folder '\t_p_' myfastint2str(gh)  '\Rad0' rad_string '\N\tijd_' tijd_string '.txt'],N_plot_5(:,:,gh));
%             else
%                 fid=fopen([save_loc folder '\t_p_' myfastint2str(gh)  '\Rad' rad_string '\N\tijd_' tijd_string '.txt'],'w');
%                 fclose(fid);
%                 csvwrite([save_loc folder '\t_p_' myfastint2str(gh)  '\Rad' rad_string '\N\tijd_' tijd_string '.txt'],N_plot_5(:,:,gh));
%             end %end if rad < 10
%             %wrting down the Velocities
%             if rad < 10
%                 fid=fopen([save_loc folder '\t_p_' myfastint2str(gh)  '\Rad0' rad_string '\V\tijd_' tijd_string '.txt'],'w');
%                 fclose(fid);
%                 csvwrite([save_loc folder '\t_p_' myfastint2str(gh)  '\Rad0' rad_string '\V\tijd_' tijd_string '.txt'],V_plot_1(:,:,gh));
%             else
%                 fid=fopen([save_loc folder '\t_p_' myfastint2str(gh)  '\Rad' rad_string '\V\tijd_' tijd_string '.txt'],'w');
%                 fclose(fid);
%                 csvwrite([save_loc folder '\t_p_' myfastint2str(gh)  '\Rad' rad_string '\V\tijd_' tijd_string '.txt'],V_plot_1(:,:,gh));
%             end %end if rad < 10
%         end %end of t_g
%         if rad < 10
%             fid=fopen([save_loc folder '\Arr\Rad0' rad_string '\N\t' tijd_string 'Ndeg.txt'],'w');
%             fclose(fid);
%             csvwrite([save_loc folder '\Arr\Rad0' rad_string '\N\t' tijd_string 'Ndeg.txt'],Ek_arr_N);
%         else
%             fid=fopen([save_loc folder '\Arr\Rad' rad_string '\N\t' tijd_string 'Ndeg.txt'],'w');
%             fclose(fid);
%             csvwrite([save_loc folder '\Arr\Rad' rad_string '\N\t' tijd_string 'Ndeg.txt'],Ek_arr_N);
%             
% 
%         end %end if rad < 10
%         clear Ek_arr_N
%     end %end of t
%     %% save plot data.
%     if sum(Total) > 0
%         
%         X_temp =[X_max_t; X_max_0; X_max_1; X_max_2; X_max_3; X_max_4];
%         P_temp =[p_plasma; p_bg];
%         T_temp =[Total_front; Total_front_deel_1; Total_front_deel_2; Total_front_deel_4 ];
%         TT_temp =[Total_front; Total_front_1; Total_front_2; Total_front_3; Total_front_4 ];
%         TTT_temp = [Total_front_M1;Total_front_M2;Total_front_M3;Total_front_M4];
%         V_temp = [V_X_max_0; V_X_max_1; V_X_max_2; V_X_max_3; V_X_max_4];
%         if rad < 10
%             %X_temp
%             fid=fopen([save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_center.txt'],'w');
%             fclose(fid);
%             csvwrite([save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_center.txt'],X_temp);
%             fid=fopen([save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_moment.txt'],'w');
%             fclose(fid);
%             csvwrite([save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_moment.txt'],P_temp);
%             fid=fopen([save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_Particles.txt'],'w');
%             fclose(fid);
%             csvwrite([save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_Particles.txt'],Np);
%             fid=fopen([save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_front.txt'],'w');
%             fclose(fid);
%             csvwrite([save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_front.txt'],T_temp);
%             fid=fopen([save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_front_inten.txt'],'w');
%             fclose(fid);
%             csvwrite([save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_front_inten.txt'],TT_temp);
%             fid=fopen([save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_front_inten.txt'],'w');
%             fclose(fid);
%             csvwrite([save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_front_inten_M.txt'],TTT_temp);
%             fid=fopen([save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_Velovity_maxima.txt'],'w');
%             fclose(fid);
%             csvwrite([save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_Velovity_maxima.txt'],V_temp);
%         else
%             fid=fopen([save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_center.txt'],'w');
%             fclose(fid);
%             csvwrite([save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_center.txt'],X_temp);
%             fid=fopen([save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_moment.txt'],'w');
%             fclose(fid);
%             csvwrite([save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_moment.txt'],P_temp);
%             fid=fopen([save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_Particles.txt'],'w');
%             fclose(fid);
%             csvwrite([save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_Particles.txt'],Np);
%             fid=fopen([save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_front.txt'],'w');
%             fclose(fid);
%             csvwrite([save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_front.txt'],T_temp);
%             fid=fopen([save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_front_inten.txt'],'w');
%             fclose(fid);
%             csvwrite([save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_front_inten.txt'],TT_temp);
%             fid=fopen([save_loc folder '\t_p_0\Rad' rad_string '\img\ana_Velovity_maxima.txt'],'w');
%             fclose(fid);
%             csvwrite([save_loc folder '\t_p_0\Rad' rad_string '\img\ana_Velovity_maxima.txt'],V_temp);
%         end %end if rad < 10
%     end
%     clear X_temp V_temp T_temp TT_temp
%     %% plot end results of rad
%     if sum(Total) > 0 && rad==Rad_delta
%         fig = figure;
%         plot(X_max_t,'xc')
%         hold on
%         plot(X_max_0,'xb')
%         plot(X_max_1,'xr')
%         plot(X_max_2,'xg')
%         hold off
%         title('Center position')
%         xlabel('timestep')
%         ylim([0 0.05])
%         legend('center plasma', '0 bots', '1 bots', '2 bots')
%         set(gca,'FontSize',16)
%         set(gca,'LineWidth',2)
%         set(get(gca,'xlabel'),'FontSize',16)
%         set(get(gca,'ylabel'),'FontSize',16)
%         set(get(gca,'title'),'FontSize',16)
%         set(findobj('Type','line'),'MarkerSize',10)
%         set(findobj('Type','line'),'LineWidth',2)
%         if rad < 10
%             saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_center.jpg'])
%             saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_center.fig'])
%         else
%             saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_center.jpg'])
%             saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_center.fig'])
%         end %end of  rad < 10
%         close
%         
%         %plot Kinetic energy
%         fig = figure;
%         semilogy(stap,Ek_plasma,'bx')
%         hold on
%         semilogy(stap,Ek_bg,'rx')
%         semilogy(stap,Ek_bg+Ek_plasma,'gx')
%         hold off
%         xlabel('tijd')
%         ylabel('Ek [J]')
%         legend('Ek_plasma', 'Ek_bg','Ek_bg+Ek_plasma')
%         title('Ek')
%         set(gca,'FontSize',16)
%         set(gca,'LineWidth',2)
%         set(get(gca,'xlabel'),'FontSize',16)
%         set(get(gca,'ylabel'),'FontSize',16)
%         set(get(gca,'title'),'FontSize',16)
%         set(findobj('Type','line'),'MarkerSize',10)
%         set(findobj('Type','line'),'LineWidth',2)
%         if rad < 10
%             saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_EK.jpg'])
%             saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_EK.fig'])
%         else
%             saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_EK.jpg'])
%             saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_EK.fig'])
%         end %end of  rad < 10
%         close
%         
%         %plot momentum
%         fig = figure;
%         semilogy(stap,p_plasma,'bx')
%         hold on
%         semilogy(stap,p_bg,'rx')
%         semilogy(stap,p_bg+p_plasma,'gx')
%         hold off
%         xlabel('tijd')
%         ylabel('momentum [m kg / s]')
%         legend('p_plasma', 'P_bg','P som')
%         title('momentum')
%         set(gca,'FontSize',16)
%         set(gca,'LineWidth',2)
%         set(get(gca,'xlabel'),'FontSize',16)
%         set(get(gca,'ylabel'),'FontSize',16)
%         set(get(gca,'title'),'FontSize',16)
%         set(findobj('Type','line'),'MarkerSize',10)
%         set(findobj('Type','line'),'LineWidth',2)
%         if rad < 10
%             saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_moment.jpg'])
%             saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_moment.fig'])
%         else
%             saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_moment.jpg'])
%             saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_moment.fig'])
%         end %end of  rad < 10
%         close
%         
%         %plot amount of particles
%         fig = figure;
%         plot(stap,Np,'bx')
%         hold on
%         %plot(stap,Nbg_st,'rx')
%         %plot(stap,Nbg_mo,'gx')
%         %plot(stap,Nbg_mo+Nbg_st,'kx')
%         hold off
%         xlabel('tijd')
%         ylabel('Amount of particles')
%         legend('N_plasma', 'Nbg_st','Nbg_mo','Nbg_som')
%         title('Particles')
%         set(gca,'FontSize',16)
%         set(gca,'LineWidth',2)
%         set(get(gca,'xlabel'),'FontSize',16)
%         set(get(gca,'ylabel'),'FontSize',16)
%         set(get(gca,'title'),'FontSize',16)
%         set(findobj('Type','line'),'MarkerSize',10)
%         set(findobj('Type','line'),'LineWidth',2)
%         if rad < 10
%             saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_Particles.jpg'])
%             saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_Particles.fig'])
%         else
%             saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_Particles.jpg'])
%             saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_Particles.fig'])
%         end %end of  rad < 10
%         close
%         
%         %plot plasma front
%         fig = figure;
%         plot(X_plot(Total_front),'cx')
%         hold on
%         plot(X_plot(round(Total_front_deel_1)),'gx')
%         plot(X_plot(round(Total_front_deel_2)),'mx')
%         plot(X_plot(round(Total_front_deel_3)),'bx')
%         plot(X_plot(round(Total_front_deel_4)),'kx')
%         hold off
%         title('Front position')
%         xlabel('timestep in 0.1us')
%         ylabel('Position')
%         xlim([0 (T_end-T_begin) / T_delta])
%         ylim([0 0.05])
%         title('Front plasma')
%         legend('Front pixel', 'Front 10% particles', 'Front 20% particles', 'Front 30% particles', 'Front 40% particles','Location','southeast')
%         set(gca,'FontSize',16)
%         set(gca,'LineWidth',2)
%         set(get(gca,'xlabel'),'FontSize',16)
%         set(get(gca,'ylabel'),'FontSize',16)
%         set(get(gca,'title'),'FontSize',16)
%         set(findobj('Type','line'),'MarkerSize',10)
%         set(findobj('Type','line'),'LineWidth',2)
%         if rad < 10
%             saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_front.jpg'])
%             saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_front.fig'])
%         else
%             saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_front.jpg'])
%             saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_front.fig'])
%         end %end of  rad < 10
%         close
%         
%         %plot plasma front
%         fig = figure;
%         plot(X_plot(Total_front+1),'cx')
%         hold on
%         plot(X_plot(Total_front_1+1),'gx')
%         plot(X_plot(Total_front_2+1),'mx')
%         plot(X_plot(Total_front_3+1),'bx')
%         plot(X_plot(Total_front_4+1),'kx')
%         hold off
%         title('Front position')
%         xlabel('timestep in 0.1us')
%         ylabel('Position')
%         xlim([0 100])
%         ylim([0 0.05])
%         title('Front plasma')
%         legend('Front pixel','Front 10% intensity', 'Front 20% intensity', 'Front 30% intensity', 'Front 40% intensity','Location','southeast')
%         set(gca,'FontSize',16)
%         set(gca,'LineWidth',2)
%         set(get(gca,'xlabel'),'FontSize',16)
%         set(get(gca,'ylabel'),'FontSize',16)
%         set(get(gca,'title'),'FontSize',16)
%         set(findobj('Type','line'),'MarkerSize',10)
%         set(findobj('Type','line'),'LineWidth',2)
%         if rad < 10
%             saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_front_inten.jpg'])
%             saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_front_inten.fig'])
%         else
%             saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_front_inten.jpg'])
%             saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_front_inten.fig'])
%         end %end of  rad < 10
%         close
%         
%         fig = figure;
%         plot(X_plot(Total_front_M+1),'cx')
%         hold on
%         plot(X_plot(Total_front_M1+1),'gx')
%         plot(X_plot(Total_front_M2+1),'mx')
%         plot(X_plot(Total_front_M3+1),'bx')
%         plot(X_plot(Total_front_M4+1),'kx')
%         hold off
%         title('Front position')
%         xlabel('timestep in 0.1us')
%         ylabel('Position')
%         xlim([0 (T_end-T_begin) / T_delta])
%         ylim([0 0.05])
%         title('Front plasma M')
%         legend('Front pixel','Front 10% intensity', 'Front 20% intensity', 'Front 30% intensity', 'Front 40% intensity','Location','southeast')
%         set(gca,'FontSize',16)
%         set(gca,'LineWidth',2)
%         set(get(gca,'xlabel'),'FontSize',16)
%         set(get(gca,'ylabel'),'FontSize',16)
%         set(get(gca,'title'),'FontSize',16)
%         set(findobj('Type','line'),'MarkerSize',10)
%         set(findobj('Type','line'),'LineWidth',2)
%         if rad < 10
%             saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_front_inten_M.jpg'])
%             saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_front_inten_M.fig'])
%         else
%             saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_front_inten_M.jpg'])
%             saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_front_inten_M.fig'])
%         end %end of  rad < 10
%         close
%         
%         fig = figure;
%         plot(X_plot(Total_front_O),'cx')
%         hold on
%         plot(X_plot(Total_front_O1),'gx')
%         plot(X_plot(Total_front_O2),'mx')
%         plot(X_plot(Total_front_O3),'bx')
%         plot(X_plot(Total_front_O4),'kx')
%         hold off
%         title('Front position')
%         xlabel('timestep in 0.1us')
%         ylabel('Position')
%         xlim([0 100])
%         ylim([0 0.05])
%         title('Front plasma O')
%         legend('Front pixel','Front 10% intensity', 'Front 20% intensity', 'Front 30% intensity', 'Front 40% intensity','Location','southeast')
%         set(gca,'FontSize',16)
%         set(gca,'LineWidth',2)
%         set(get(gca,'xlabel'),'FontSize',16)
%         set(get(gca,'ylabel'),'FontSize',16)
%         set(get(gca,'title'),'FontSize',16)
%         set(findobj('Type','line'),'MarkerSize',10)
%         set(findobj('Type','line'),'LineWidth',2)
%         if rad < 10
%             saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_front_inten_O.jpg'])
%             saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_front_inten_O.fig'])
%         else
%             saveas(fig,[save_loc folder '\t_p_0'\Rad' rad_string '\img\ana_peak_front_inten_O.jpg'])
%             saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_front_inten_O.fig'])
%         end %end of  rad < 10
%         close
%         
%         
%         fig = figure;
%         plot(X_plot(Total_front_O_2),'cx')
%         hold on
%         plot(X_plot(Total_front_M_2),'gx')
%         hold off
%         title('Front position')
%         xlabel('timestep in 0.1us')
%         ylabel('Position')
%         xlim([0 100])
%         ylim([0 0.05])
%         title('Front plasma MO')
%         legend('Front O','Front M')
%         set(gca,'FontSize',16)
%         set(gca,'LineWidth',2)
%         set(get(gca,'xlabel'),'FontSize',16)
%         set(get(gca,'ylabel'),'FontSize',16)
%         set(get(gca,'title'),'FontSize',16)
%         set(findobj('Type','line'),'MarkerSize',10)
%         set(findobj('Type','line'),'LineWidth',2)
%         if rad < 10
%             saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_front_inten_MO.jpg'])
%             saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_peak_front_inten_MO.fig'])
%         else
%             saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_front_inten_MO.jpg'])
%             saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\ana_peak_front_inten_MO.fig'])
%         end %end of  rad < 10
%         close
%         
%         fig = figure;
%         plot(X_max_t,'cx')
%         hold on
%         plot(X_max_0,'bx')
%         plot(X_max_1,'rx')
%         plot(X_max_2,'gx')
%         plot(X_max_3,'mx')
%         plot(X_max_4,'yx')
%         hold off
%         title('Positie Maxima')
%         ylim([0 0.05])
%         xlabel('timestep in 0.1us')
%         ylabel('Distance [m]')
%         legend('Plasma','0bots', '1bots','2bots','3bots')
%         set(gca,'FontSize',16)
%         set(gca,'LineWidth',2)
%         set(get(gca,'xlabel'),'FontSize',16)
%         set(get(gca,'ylabel'),'FontSize',16)
%         set(get(gca,'title'),'FontSize',16)
%         set(findobj('Type','line'),'MarkerSize',10)
%         set(findobj('Type','line'),'LineWidth',2)
%         if rad < 10
%             saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_maxima_position.jpg'])
%             saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_maxima_position.fig'])
%         else
%             saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\ana_maxima_position.jpg'])
%             saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\ana_maxima_position.fig'])
%         end %end of  rad < 10
%         close
%         
%         fig = figure;
%         plot(V_X_max_0,'bx')
%         hold on
%         plot(V_X_max_1,'rx')
%         plot(V_X_max_2,'gx')
%         plot(V_X_max_3,'mx')
%         plot(V_X_max_4,'yx')
%         hold off
%         title('Velocity maxima')
%         xlabel('timestep in 0.1us')
%         ylabel('Velocity in [m/s]')
%         legend('0bots', '1bots','2bots','3bots','4bots')
%         set(gca,'FontSize',16)
%         set(gca,'LineWidth',2)
%         set(get(gca,'xlabel'),'FontSize',16)
%         set(get(gca,'ylabel'),'FontSize',16)
%         set(get(gca,'title'),'FontSize',16)
%         set(findobj('Type','line'),'MarkerSize',10)
%         set(findobj('Type','line'),'LineWidth',2)
%         if rad < 10
%             saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_Velovity_maxima.jpg'])
%             saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_Velovity_maxima.fig'])
%         else
%             saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\ana_Velovity_maxima.jpg'])
%             saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\ana_Velovity_maxima.fig'])
%         end %end of  rad < 10
%         close
%         
%         %EK_arr
%         fig = figure;
%         %plot(sum(Ek_arr_V.^2*m_c),'bx')
%         hold on
%         % plot(Ek_arr(1,t),'bx')
%         % plot(Ek_arr(2,t),'bx')
%         hold off
%         title('Kinetic energy arriving species')
%         xlabel('timestep in 0.1us')
%         ylabel('Ek in [J]')
%         legend('Total_Ek', 'Ek 0bots','Ek 1bots','3bots','4bots')
%         set(gca,'FontSize',16)
%         set(gca,'LineWidth',2)
%         set(get(gca,'xlabel'),'FontSize',16)
%         set(get(gca,'ylabel'),'FontSize',16)
%         set(get(gca,'title'),'FontSize',16)
%         set(findobj('Type','line'),'MarkerSize',10)
%         set(findobj('Type','line'),'LineWidth',2)
%         if rad < 10
%             saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_Ek_arr_t_maxima.jpg'])
%             saveas(fig,[save_loc folder '\t_p_0\Rad0' rad_string '\img\ana_Ek_arr_t_maxima.fig'])
%         else
%             saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\ana_Velovity_maxima.jpg'])
%             saveas(fig,[save_loc folder '\t_p_0\Rad' rad_string '\img\ana_Velovity_maxima.fig'])
%         end %end of  rad < 10
%         close
    end %if total>1
    
%     %% write data into files
%     %write sum values for collisions
%     if rad < 10
%         fid=fopen([save_loc folder '\t_p_0\Rad0' rad_string 'deg.txt'],'w');
%         fclose(fid);
%         csvwrite([save_loc folder '\t_p_0\Rad0' rad_string 'deg.txt'],Totaal_t);
%     else
%         fid=fopen([save_loc folder '\t_p_0\Rad' rad_string 'deg.txt'],'w');
%         fclose(fid);
%         csvwrite([save_loc folder '\t_p_0\Rad' rad_string 'deg.txt'],Totaal_t);
%     end %end if rad < 10
    
end %end of rad
% All data points have been generated
%% Solve for the 3D plume
%--------------------------------------------------------------------------


% close all
% clear M M_1 Mbg N_plot N_plot_5 V_plot V_plot_1 Totaal_t Total_front Total_front_1 Total_front_2 Total_front_3 n_nbg
% clear N_loc_2 N_max_2 k Nbg_plot Total Nbg_plot_1 Nbg_Vplot_1 Nbg_Vplot_3 S_p S_exp S_mo s_size Sbg_mo Sbg_nb Sbg_st Np nor_V_temp ssx ssx_1 ssx_2
% clear V_plot_1 N_plot_2 a N_plot V_plot N_plot_5 Nbg_plot_2 Nbg_Vplot_3 X_plot nor_V O2_bg p p_bg P_dep oo stap t  t_oo ttg t_o t_pp
% clear theta k_plot Np O2_bg Nbg_plot_som Vol_T v_size x_names high high_T Ev K Ek_V Ek_Vc gh k_numb lost Max_x_size Nbg_mo Nbg_plot
% clear X_n gg gh bx cos_theta dis_iPoint dis_vPlotMax X_max_0 X_max_1 X_max_2 X_max_3 X_max_4 X_max_t X_n i_nb_p Ek_arr_N Ek_bg Ek_loc fig
% clear Vbg Vbg_new Vbg_min_st dis_E_x comment i z xx ttijd ttf temp ssx N_loc_1 N_loc_0 N_loc_3 V_X_max_0 init_N_dist init_V_dist K X V
% clear S_exp S_mo S_new_nbots S_st Total_front_4 Total_front_deel_3 Total_front_deel_2 Total_front_deel_1 Nbg_st nn nnb N_diff_bg N_diff_p N_loc
% clear V P_dep Ek_N N p P P_temp V_plot_2 N_loc_1 N_loc_2 N_loc_3 N_loc_4 N_loc_t N_loc_0 N_max_0 N_max_1 N_max_2 N_max_3 N_max_t p_bg P_local S_exp S_mo
% clear Total_front_1 Total_front Total_front_3 Total_front_2 Total_front_4 Total_front_deel_3 Total_front_deel_1 Total_front_deel_2 V_X_max_0
% clear V_X_max_1 V_X_max_2 V_X_max_3 V_X_max_4 V X_n_b_mo X_n_b_st X_n_nb fid kb loc lost a b R_kans_botsen_st_exp R_kans_botsen_mo_exp Total_front_deel_4
% clear theta A_coeff Vol_T t S_newnbotsbg R_kans_bots_mo R_kans_bots_st rad p_plasma V_sort_loc V_sort_val V_sort_loc V_sort_val Ek_plasma
% clear dis_v_dist_x dis_N_dist_x N_atom_temp X_nbg_b_st_a x_new_bots_mo v_gem VOLbg Vbg_min_mo size_EK size_norv size_V_dist N_atoom_rad Nbg Nbg_exp
% clear Xbg_n_minb Xbg_n_minb_curr Xbg_n_b_st Xbg_n_b_mo Xbg_n_b_mo T_sum Nbg_exp_N Nbg_exp_V Nbg_size mm sm x_new_nbotsbg x_newbg_bots_mo x_newbg_bots_st
% clear x_newbg_nbots R_kans_nbots ip jj m Kbg K_local i_b_k i_b_p einde_N fitplotN fitplotV MO2 N_atoom N_atoom_totaal fitplotX firplot Y
% clear sizeN5 testN testV V_plot_7 Total_O Total_M Total_O2 dd x vv smooth_fraction Total_front_M Total_front_M1 Total_front_M2 Total_front_M3 Total_front_M4
% clear Total_front_M_2 Total_front_O Total_front_O1 Total_front_O2 Total_front_O3 Total_front_O4 Total_front_O_2 Total_front_M1 Total_front_M Total_front_M2
% clear Total_front_M3 Total_front_M4 Spot_size Temperature_BG tijd_string value_X_n_b_mo_a value_X_nbg_b_mo_a value_X_n_b_st value_X_nbg_b_st_a Vbg_exp
% clear names N_temp_left N_loc_b N_loc_nb loc_particles loc_Velocity make_1d_plots X_plot_5_lim X_str_end1 X_str_end2 Xbg_n_nb N_plot_7 Min_angle_particles
% clear o nor save_loc Total_Ek_theory Total_p_theory V_min_mo V_min_st V_stepsize i_b_v i_nb_k i_nb_v limbots loc_b loc_nb n ans V_fit Oxidize
% clear Opp_plasma Begin_N bots density_bg Density_target dis_NPoints E_O_high fitplotY i_b_size i_nb_size yl ym x_new_nbots x_new_bots_st P_local_temp
% clear part E_O_low V_temp V_maximaal V_local_dist theta_i E_O_low Eb_UC Rad_stop show_BG_inplot rad_string p_plot m_c m_o N_local_dist N_local_dist_min
% 
% 
% if dried == 1
%     abb = Plot_3d_full(plek ,X_delta, Rad_delta,X_end,(2),T_end,N_atoom_totaal_N,tijden,X_TS)
% end %end of dried==1
% if Arr_sum == 1
%     abb  = Plot_Arr_full(plek, Rad_delta, t_p, T_end, m_p, V_minimaal, folder, T_delta, X_begin, X_end, X_delta)
% end
% profile off
