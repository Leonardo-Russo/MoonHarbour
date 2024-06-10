%% Ephemeris Visualizer - Leonardo Russo

close all
clear
clc
                                                                                
addpath('Library/')
addpath('Data/')
addpath('Data/Planets/')
addpath('Data/Materials/')
addpath('Data/Ephemeris/')

% Define Options
global opt
opt = struct('name', "Options");
opt.saveplots = false;
opt.create_animation = false;
opt.show_progress = false;
opt.compute_target = true;
opt.compute_direct_approach = true;
opt.include_actuation = false;
opt.additional_plots = false;
opt.showgui = false;
opt.N = 1000;                   % n° of points for the Interpolation
opt.RelTolODE = 1e-7;           % options for ode()
opt.AbsTolODE = 1e-6;

OptionsODE = odeset('RelTol', opt.RelTolODE, 'AbsTol', opt.AbsTolODE);

% Define Global Variables
global DU TU MU Rm muM pbar log

% Define Constants
DU = 1738;                                              % km
TU = sqrt(1 / 4902.7779 * DU^3);                        % s
Rm = 1738.1 / DU;                                       % DU
muM = 4902.7779 * TU^2 / DU^3;                          % DU^3/TU^2
muE = 398600.4418 * TU^2 / DU^3;                        % DU^3/TU^2
muS = 132712440018 * TU^2 / DU^3;                       % DU^3/TU^2
deltaE = deg2rad(23.45);                                % rad
psiM = deg2rad(-81.7 + 360/18.6 * (5 + 4.5/12));        % rad
deltaM = deg2rad(1.5);                                  % rad
Dsol = 86400;                                           % s
sec2hrs = 1/3600;                                       % hrs
MU = 1/DU^2;

% Define Time Domain
date0 = datetime('2025-05-23 9:56:10');
datef = datetime('2025-05-30 9:56:10');    % @ periselenium

% Define the Chaser Initial Conditions ~ norm(rhodot_LVLH) = 1 m/s -> condition applied for the finetuning of the gain parameters such as time deltas
rho0_LVLH = [1.5, 0, 0]'/DU;                        % km
rhodot0_LVLH = [-1e-6, -1e-3, -1e-3]'/DU*TU;        % km/s
RHO0_LVLH = [rho0_LVLH; rhodot0_LVLH];

% Define Direct Approach Conditions
final_dist = 5e-3/DU;       % 5 m
final_vel = -5e-6/DU*TU;    % -5 mm/s
rhof_LVLH_DA = final_dist * RHO0_LVLH(1:3)/norm(RHO0_LVLH(1:3));
rhodotf_LVLH_DA = final_vel * RHO0_LVLH(1:3)/norm(RHO0_LVLH(1:3));
RHOf_LVLH_DA = [rhof_LVLH_DA; rhodotf_LVLH_DA];


%% Show Reference Trajectory

everything = 1;

% Load Data from Ephemeris - states are provided in ECI
load('Ephemeris.mat', 'stateDSG', 'stateMoon', 'stateSun', 'time');

if everything

    Npoints = size(stateDSG, 1);
    
    % Retrieve Data from Ephemeris
    stateDSG_ECI = stateDSG;
    stateMoon_ECI = stateMoon;
    stateSun_ECI = stateSun;
    stateEarth_ECI = zeros(Npoints, 6);

else

    % Look for Initial and Final times
    t0_sharp = juliandate(date0)*Dsol;
    tf_sharp = juliandate(datef)*Dsol;
    
    t0 = interp1(time, time, t0_sharp, 'nearest');
    t0_idx = find(time == t0);
    
    tf = interp1(time, time, tf_sharp, 'nearest');
    tf_idx = find(time == tf); 

    t0_idx = 1;
    tf_idx = 60000;
    
    % Retrieve Data from Ephemeris
    stateDSG_ECI = stateDSG(t0_idx : tf_idx, :);
    stateMoon_ECI = stateMoon(t0_idx : tf_idx, :);
    stateSun_ECI = stateSun(t0_idx : tf_idx, :);
    stateEarth_ECI = zeros(length(t0_idx : tf_idx), 6);
    time = time(t0_idx : tf_idx);

    Npoints = length(time);

end

% Conversion into Canonical Units
stateDSG_ECI(:, 1:3) = stateDSG_ECI(:, 1:3) / DU;
stateDSG_ECI(:, 4:6) = stateDSG_ECI(:, 4:6) / DU * TU;
stateMoon_ECI(:, 1:3) = stateMoon_ECI(:, 1:3) / DU;
stateMoon_ECI(:, 4:6) = stateMoon_ECI(:, 4:6) / DU * TU;
stateSun_ECI(:, 1:3) = stateSun_ECI(:, 1:3) / DU;
stateSun_ECI(:, 4:6) = stateSun_ECI(:, 4:6) / DU * TU;
stateEarth_ECI(:, 1:3) = stateEarth_ECI(:, 1:3) / DU;
stateEarth_ECI(:, 4:6) = stateEarth_ECI(:, 4:6) / DU * TU;

% Initialize Local Variables
stateEarth_MCI = zeros(Npoints, 6);
stateDSG_MCI = zeros(Npoints, 6);
stateSun_MCI = zeros(Npoints, 6);

time = time/TU;
t0 = time(1);
tf = time(end);

for i = 1 : Npoints

    % Retrieve ECI Position and Velocity Vectors
    rDSG_ECI = stateDSG_ECI(i, 1:3)';
    vDSG_ECI = stateDSG_ECI(i, 4:6)';

    rM_ECI = stateMoon_ECI(i, 1:3)';
    vM_ECI = stateMoon_ECI(i, 4:6)';

    rS_ECI = stateSun_ECI(i, 1:3)';
    vS_ECI = stateSun_ECI(i, 4:6)';

    rE_ECI = stateEarth_ECI(i, 1:3)';
    vE_ECI = stateEarth_ECI(i, 4:6)';

    % Define ECI2MCI Rotation Matrices
    R1deltaE = R1(-deltaE);
    R3PsiM = R3(psiM);
    R1deltaM = R1(deltaM);

    % Perform ECI2MCI Conversion
    rDSG_MCI = ((rDSG_ECI - rM_ECI)' * R1deltaE * R3PsiM' * R1deltaM')';
    vDSG_MCI = ((vDSG_ECI - vM_ECI)' * R1deltaE * R3PsiM' * R1deltaM')';

    rE_MCI = ((rE_ECI - rM_ECI)' * R1deltaE * R3PsiM' * R1deltaM')';
    vE_MCI = ((vE_ECI - vM_ECI)' * R1deltaE * R3PsiM' * R1deltaM')';

    rS_MCI = ((rS_ECI - rM_ECI)' * R1deltaE * R3PsiM' * R1deltaM')';
    vS_MCI = ((vS_ECI - vM_ECI)' * R1deltaE * R3PsiM' * R1deltaM')';

    % Compute Initial State
    if i == 1
        X0t_MCI = [rDSG_MCI; vDSG_MCI];       % DSG State in MCI
        COE0t = rvPCI2COE(X0t_MCI', muM);
        MEE0t = COE2MEE(COE0t)';
    end


    % Assign values to the Matrices
    stateEarth_MCI(i, :) = [rE_MCI', vE_MCI'];
    stateSun_MCI(i, :) = [rS_MCI', vS_MCI'];
    stateDSG_MCI(i, :) = [rDSG_MCI', vDSG_MCI'];

end

% Perform the Interpolation
EarthPPsMCI = get_statePP(time, stateEarth_MCI);
SunPPsMCI = get_statePP(time, stateSun_MCI);
% DSGPPsMCI = get_statePP(time, stateDSG_MCI);
MoonPPsECI = get_statePP(time, stateMoon_ECI);

fig_MCI = figure('Name', 'Traj');
DrawTrajMCI3D(stateDSG_MCI*DU, '#03adfc');
savefig(fig_MCI, 'reference_trajectory.fig');
print(fig_MCI, 'reference_trajectory.png', '-dpng', '-r600');          % 300 DPI
view([0, 90]);
print(fig_MCI, 'reference_trajectory_top.png', '-dpng', '-r600');          % 300 DPI

% Define the timespan for the propagation
tspan_ref = time;

% Perform the Propagation of the Target Trajectory
[~, MEEt_ref] = ode113(@(t, MEE) DynamicalModelMEE(t, MEE, EarthPPsMCI, SunPPsMCI, ...
    muE, muS, MoonPPsECI, deltaE, psiM, deltaM, t0, tf), tspan_ref, MEE0t, OptionsODE);

% Conversion from MEE to COE
COEt_ref = MEE2COE(MEEt_ref);

% Conversion from COE to MCI
Xt_MCI_ref = COE2rvPCI(COEt_ref, muM);

% % Interpolate the Angular Velocity of LVLH wrt MCI
% [omegaPPsLVLH, omegadotPPsLVLH] = TargetHandler(Xt_MCI_ref, COEt_ref, MEEt_ref, tspan_ref, ...
%     EarthPPsMCI, SunPPsMCI, MoonPPsECI, deltaE, psiM, deltaM, muE, muS);


% Transform into Synodic Frame
stateDSG_SYN = MCI2SYN(stateDSG_MCI, time, MoonPPsECI, deltaE, psiM, deltaM, 1, stateMoon_ECI);

% Plot the trajectory in the synodic frame
fig_syn = figure('Name', 'Trajectory in Synodic Frame');
DrawTrajMCI3D(stateDSG_SYN(:, 1:3)*DU, '#03adfc');
savefig(fig_syn, 'reference_trajectory_SYN.fig');
print(fig_syn, 'reference_trajectory_SYN.png', '-dpng', '-r600');

view([0, 90]);
print(fig_syn, 'reference_trajectory_SYN_top.png', '-dpng', '-r600');



% % Compute COE for DSG
% COE_DSG = rvPCI2COE(stateDSG_MCI, muM);
% 
% % Confront the Classical Orbital Elements
% a_fig = figure('name', 'a Comparison');
% plot((time-t0)*TU/Dsol, COE_DSG(:, 1)*DU, 'LineStyle','-', 'LineWidth', 1.5, 'Color', '#03adfc')
% hold on
% plot((time-t0)*TU/Dsol, COEt_ref(:, 1)*DU, 'LineStyle','-', 'LineWidth', 1.5, 'Color', '#fc9803')
% grid on
% xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 12)
% ylabel('$a \ [km]$', 'Interpreter','latex', 'FontSize', 12)
% legend('Predicted', 'Reference', 'location', 'best')
% savefig(a_fig, 'coe_stuff/a_comp.fig');
% print(a_fig, 'coe_stuff/a_comp.png', '-dpng', '-r600');
% 
% e_fig = figure('name', 'e Comparison');
% plot((time-t0)*TU/Dsol, COE_DSG(:, 2), 'LineStyle','-', 'LineWidth', 1.5, 'Color', '#03adfc')
% hold on
% plot((time-t0)*TU/Dsol, COEt_ref(:, 2), 'LineStyle','-', 'LineWidth', 1.5, 'Color', '#fc9803')
% grid on
% xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 12)
% ylabel('$e$', 'Interpreter','latex', 'FontSize', 12)
% legend('Predicted', 'Reference', 'location', 'best')
% savefig(e_fig, 'coe_stuff/e_comp.fig');
% print(e_fig, 'coe_stuff/e_comp.png', '-dpng', '-r600');
% 
% i_fig = figure('name', 'i Comparison');
% plot((time-t0)*TU/Dsol, rad2deg(COE_DSG(:, 3)), 'LineStyle','-', 'LineWidth', 1.5, 'Color', '#03adfc')
% hold on
% plot((time-t0)*TU/Dsol, rad2deg(COEt_ref(:, 3)), 'LineStyle','-', 'LineWidth', 1.5, 'Color', '#fc9803')
% grid on
% xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 12)
% ylabel('$i \ [deg]$', 'Interpreter','latex', 'FontSize', 12)
% legend('Predicted', 'Reference', 'location', 'best')
% savefig(i_fig, 'coe_stuff/i_comp.fig');
% print(i_fig, 'coe_stuff/i_comp.png', '-dpng', '-r600');
% 
% Omega_fig = figure('name', 'Omega Comparison');
% plot((time-t0)*TU/Dsol, rad2deg(COE_DSG(:, 4)), 'LineStyle','-', 'LineWidth', 1.5, 'Color', '#03adfc')
% hold on
% plot((time-t0)*TU/Dsol, rad2deg(COEt_ref(:, 4)), 'LineStyle','-', 'LineWidth', 1.5, 'Color', '#fc9803')
% grid on
% xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 12)
% ylabel('$\Omega \ [deg]$', 'Interpreter','latex', 'FontSize', 12)
% legend('Predicted', 'Reference', 'location', 'best')
% savefig(Omega_fig, 'coe_stuff/Omega_comp.fig');
% print(Omega_fig, 'coe_stuff/Omega_comp.png', '-dpng', '-r600');
% 
% omega1_fig = figure('name', 'omega Comparison');
% plot((time-t0)*TU/Dsol, rad2deg(COE_DSG(:, 5)), 'LineStyle','-', 'LineWidth', 1.5, 'Color', '#03adfc')
% hold on
% plot((time-t0)*TU/Dsol, rad2deg(COEt_ref(:, 5)), 'LineStyle','-', 'LineWidth', 1.5, 'Color', '#fc9803')
% grid on
% xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 12)
% ylabel('$\omega \ [deg]$', 'Interpreter','latex', 'FontSize', 12)
% legend('Predicted', 'Reference', 'location', 'best')
% savefig(omega1_fig, 'coe_stuff/omegapiccolo_comp.fig');
% print(omega1_fig, 'coe_stuff/omegapiccolo_comp.png', '-dpng', '-r600');
% 
% theta_fig = figure('name', 'theta Comparison');
% plot((time-t0)*TU/Dsol, rad2deg(COE_DSG(:, 6)), 'LineStyle','-', 'LineWidth', 1.5, 'Color', '#03adfc')
% hold on
% plot((time-t0)*TU/Dsol, rad2deg(COEt_ref(:, 6)), 'LineStyle','-', 'LineWidth', 1.5, 'Color', '#fc9803')
% grid on
% xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 12)
% ylabel('$\theta^{*} \ [deg]$', 'Interpreter','latex', 'FontSize', 12)
% legend('Predicted', 'Reference', 'location', 'best')
% savefig(theta_fig, 'coe_stuff/theta_comp.fig');
% print(theta_fig, 'coe_stuff/theta_comp.png', '-dpng', '-r600');


