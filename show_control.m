%% Moon Harbour Project - Leonardo Russo

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

sampling_time = 30;                     % seconds
verbose = true;
misalignment_type = "null";
state_perturbation_flag = false;
engine_failure_flag = false;
workspace_path = 'Data/Utils/show_control.mat';

OptionsODE = odeset('RelTol', opt.RelTolODE, 'AbsTol', opt.AbsTolODE);

% Define Thrust Misalignment
misalignments = [];
misalignment0 = define_misalignment_error(misalignment_type);
misalignments = [misalignments; misalignment0];


tic
%% Hyperparameters and Settings

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
datef = datetime('2025-05-23 21:56:10');    % @ periselenium
total_time = datef - date0;

% Define Control Limit
g0 = 9.80665;           % m/s^2
u_limit = 5e-5*g0;      % m/s^2
u_lim = u_limit * 1e-3 * TU^2 / DU;

% Define the Chaser Initial Conditions ~ norm(rhodot_LVLH) = 1 m/s -> condition applied for the finetuning of the gain parameters
RHO0_LVLH = [1.5, 0, 0, -1e-6, -1e-3, -1e-3]';              % km, km/s
RHO0_LVLH = [RHO0_LVLH(1:3)/DU; RHO0_LVLH(4:6)/DU*TU];      % adim

% Define Direct Approach Conditions
rhof_LVLH_DA = 5e-3/DU * RHO0_LVLH(1:3)/norm(RHO0_LVLH(1:3));
rhodotf_LVLH_DA = -5e-5/DU*TU * RHO0_LVLH(1:3)/norm(RHO0_LVLH(1:3));
RHOf_LVLH_DA = [rhof_LVLH_DA; rhodotf_LVLH_DA];


%% Propagate Reference Target Trajectory

% Interpolate the Ephemeris and Retrieve Target's Initial State
[X0t_MCI, COE0t, MEE0t, EarthPPsMCI, DSGPPsMCI, SunPPsMCI, MoonPPsECI, time, t0, tf, opt.N] = ...
    EphemerisHandlerExp(deltaE, psiM, deltaM, opt.N, date0, datef);

% Define the timespan for the propagation
tspan_ref = linspace(t0, tf, opt.N);
dt_ref = tspan_ref(2) - tspan_ref(1);

% Perform the Propagation of the Target Trajectory
if opt.show_progress
    pbar = waitbar(0, 'Performing the Target Trajectory Propagation');
end
[~, MEEt_ref] = ode113(@(t, MEE) DynamicalModelMEE(t, MEE, EarthPPsMCI, SunPPsMCI, ...
    muE, muS, MoonPPsECI, deltaE, psiM, deltaM, t0, tf), tspan_ref, MEE0t, OptionsODE);
if opt.show_progress
    close(pbar);
end

% Conversion from MEE to COE
COEt_ref = MEE2COE(MEEt_ref);

% Conversion from COE to MCI
Xt_MCI_ref = COE2rvPCI(COEt_ref, muM);

% Interpolate the Angular Velocity of LVLH wrt MCI
[~, omegadotPPsLVLH] = TargetHandler(Xt_MCI_ref, COEt_ref, MEEt_ref, tspan_ref, ...
    EarthPPsMCI, SunPPsMCI, MoonPPsECI, deltaE, psiM, deltaM, muE, muS);
    

%% Direct Approach: Chaser Backwards Propagation from Final Conditions

% Combine the Target and Chaser States into TC ~ [Target; Chaser] State
TC0 = [MEE0t; RHO0_LVLH];

% Define Natural Drift Propagation Settings
optODE_back = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'Events', @(t, Y) driftstop_back(t, Y));

% Define the Backwards Drift tspan
tspan_back_DA = linspace(tf, t0 + (tf-t0)/2, opt.N);

% Define Final Conditions
MEEft = MEEt_ref(end, :)';
TCf_backdrift_DA = [MEEft; RHOf_LVLH_DA];

% Perform the Backpropagation of Target and Chaser Trajectories
[tspan_back_DA, TC_backdrift_DA, t_bwdstop_DA, ~, ~] = ode113(@(t, TC) NaturalRelativeMotion(t, TC, EarthPPsMCI, SunPPsMCI, ...
    muE, muS, MoonPPsECI, deltaE, psiM, deltaM, t0, tf, omegadotPPsLVLH, 0), tspan_back_DA, TCf_backdrift_DA, optODE_back);

% Check Successful Back Drift
if isempty(t_bwdstop_DA)
    warning('Could not complete Backwards Natural Drift from Direct Approach Conditions.');
end

% Retrieve the States and epoch for Final Drift
TC0_backdrift_DA = TC_backdrift_DA(end, :)';
t0_backdrift_DA = tspan_back_DA(end);


%% Direct Approach: Generate the Reference Chaser Trajectory for Hybrid Predictive Control

% Set the Via Points and Interpolate the Reference Trajectory
[RHOdPPsLVLH_DA, viapoints_DA, t_viapoints_DA] = ReferenceTrajectory(TC0, TC0_backdrift_DA, t0, t0_backdrift_DA, [0, 1]);

% Control Settings
prediction_interval_DA = 120;                  % s - one prediction every 2 minutes
prediction_dt_DA = prediction_interval_DA/TU;     % adim - prediction interval adimensionalized

% prediction_maxlength_DA = fix(seconds(total_time)/prediction_interval_DA);        % max length of the predictive controls
prediction_maxlength_DA = 2;
check_times_DA = zeros(prediction_maxlength_DA, 1);           % gets all the times at which one should check in continue to saturate by future prediction

for k = 1 : prediction_maxlength_DA
    check_times_DA(k) = t0 + (k-1) * prediction_dt_DA;        % compute the check_times vector
end

prediction_delta_DA = (30*60)/TU;      % predictive propagation of 30 min forward
N_inner_DA = round(prediction_delta_DA/dt_ref) - 1;       % n° of points in the inner propagation


%% Direct Approach: Perform Chaser Rendezvous Manoeuvre

% Define Initial Conditions
x70 = 1;                % initial mass ratio
TCC0 = [TC0; x70];      % initial TargetControlledChaser State
event_odefun = 1;       % 1 means it'll stop at saturation

% [tspan_ctrl_DA, TCC_ctrl_DA] = odeHamHPC(@(t, TCC) AdaptableFeedbackControl(t, TCC, EarthPPsMCI, SunPPsMCI, muM, ...
%     muE, muS, MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, RHOdPPsLVLH_DA, DU, TU),...
%     [t0, t0_backdrift_DA], TCC0, opt.N);

[tspan_ctrl_DA, TCC_ctrl_DA] = odeHamHPC(@(t, TCC) HybridPredictiveControl(t, TCC, EarthPPsMCI, SunPPsMCI, muM, ...
    muE, muS, MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, tf, RHOdPPsLVLH_DA, u_lim, DU, TU, check_times_DA, prediction_delta_DA, N_inner_DA, misalignment0),...
    [t0, t0_backdrift_DA], TCC0, opt.N, @noevent, event_odefun);

% [tspan_ctrl_DA, TCC_ctrl_DA] = odeHamHPC(@(t, TCC) FixedSaturationControl(t, TCC, EarthPPsMCI, SunPPsMCI, muM, ...
%     muE, muS, MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, tf, RHOdPPsLVLH_DA, u_lim, DU, TU, check_times_DA, prediction_delta_DA, N_inner_DA, misalignment0),...
%     [t0, t0_backdrift_DA], TCC0, opt.N, @noevent, event_odefun);

% Define stopSaturationTime
stopSaturationTime = tspan_ctrl_DA(end);

% Retrive Final Kp
[~, ~, ~, ~, ~, ~, ~, ~, ~, ~, kp, ~] = HybridPredictiveControl(tspan_ctrl_DA(end), TCC_ctrl_DA(end, :), EarthPPsMCI, SunPPsMCI, muM, ...
    muE, muS, MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, tf, RHOdPPsLVLH_DA, u_lim, DU, TU, check_times_DA, prediction_delta_DA, N_inner_DA, misalignment0);

% Perform Natural Feedback Control
[tspan_temp, TCC_temp] = odeHamHPC(@(t, TCC) NaturalFeedbackControl(t, TCC, EarthPPsMCI, SunPPsMCI, muM, ...
    muE, muS, MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, tf, RHOdPPsLVLH_DA, kp, u_lim, DU, TU, misalignment0, opt.show_progress, 0), ...
    [tspan_ctrl_DA(end), t0_backdrift_DA], TCC_ctrl_DA(end, :), opt.N, @noevent);

% Retrieve Final State Values
tspan_ctrl_DA = [tspan_ctrl_DA; tspan_temp];
TCC_ctrl_DA = [TCC_ctrl_DA; TCC_temp];

% Check if Terminal Conditions are reached
[~, terminal_flag] = is_terminal_distance(tspan_ctrl_DA(end), TCC_ctrl_DA(end, :));
if ~terminal_flag
    error("Didn't reach the end...")
end

M_ctrl_DA = length(tspan_ctrl_DA);

u_temp = zeros(M_ctrl_DA, 3);
u_norm_temp = zeros(M_ctrl_DA, 1);
f_norm_temp = zeros(M_ctrl_DA, 1);
kp_temp = zeros(M_ctrl_DA, 1);

for k = 1 : M_ctrl_DA

    % [~, ~, ~, ~, u_temp(k, :), ~, ~, ~, f_norm_temp(k), kp_temp(k)] = AdaptableFeedbackControl(tspan_ctrl_DA(k), TCC_ctrl_DA(k, :), EarthPPsMCI, SunPPsMCI, muM, ...
    % muE, muS, MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, RHOdPPsLVLH_DA, DU, TU);
    
    [~, ~, ~, ~, u_temp(k, :), ~, ~, ~, ~, f_norm_temp(k), kp_temp(k), ~] = ...
        HybridPredictiveControlPostProcessing(tspan_ctrl_DA(k), TCC_ctrl_DA(k, :), EarthPPsMCI, SunPPsMCI, muM, ...
        muE, muS, MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, tf, RHOdPPsLVLH_DA, DU, TU, stopSaturationTime, misalignment0);
    
    u_norm_temp(k) = norm(u_temp(k, :));

end

% Control Norm
fig = figure('name', 'Control Thrust');
p1 = plot((tspan_ctrl_DA-t0)*TU*sec2hrs, u_norm_temp*1000*DU/TU^2, 'Color', '#4195e8', 'LineWidth', 1.5);
hold on
p2 = plot((tspan_ctrl_DA-t0)*TU*sec2hrs, u_limit*ones(length(tspan_ctrl_DA), 1), 'r--', 'LineWidth', 1.2);
p3 = plot((tspan_ctrl_DA-t0)*TU*sec2hrs, f_norm_temp*1000*DU/TU^2, 'Color', '#93faad', 'LineWidth', 1.5);
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$[m/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
title('Control Norm')
grid on
legend([p1, p2, p3], '$|u|$', '$u_{max}$', '$f$','Location', 'best', 'Fontsize', 12, 'Interpreter', 'latex');

% savefig(fig, 'natural_control_norm.fig');
% print(fig, 'natural_control_norm.png', '-dpng', '-r300');          % 300 DPI

savefig(fig, 'hp_control_norm2.fig');
print(fig, 'hp_control_norm2.png', '-dpng', '-r300');          % 300 DPI

