%% Moon Harbour Project - Leonardo Russo

close all
clear
clc
                                                                                
addpath('Library/')
addpath('Library/Dario/')
addpath('Data/')
addpath('Data/Planets/')
addpath('Data/Materials/')
addpath('Data/temp/')


% Introduce Options Structure
opt = struct('name', "Progator Options");
opt.saveplots = false;
opt.create_animation = false;

% Define options for ode113()
opt.RelTolODE = 1e-7;
opt.AbsTolODE = 1e-6;
OptionsODE = odeset('RelTol', opt.RelTolODE, 'AbsTol', opt.AbsTolODE);


%% Hyperparameters and Settings

% Define Global Variables
global DU TU Rm muM pbar log

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

% Define Time Domain
date0 = datetime('2025-05-23 9:56:10');
datef = datetime('2025-05-23 21:56:10');    % @ periselenium
total_time = datef - date0;

% Define the Chaser Initial Conditions ~ norm(rhodot_LVLH) = 1 m/s
% This was the condition applied for the finetuning of the gain parameters
RHO0_LVLH = [500e-3, 0, 0, 0, -1e-4, -1e-4]';                      % km, km/s
RHO0_LVLH = [RHO0_LVLH(1:3)/DU; RHO0_LVLH(4:6)/DU*TU];      % adim

% Define Desired Conditions for Docking
RHOf_LVLH = [-5e-3, 0, 0, 1e-5, 0, 0]';                 % km, km/s
RHOf_LVLH = [RHOf_LVLH(1:3)/DU; RHOf_LVLH(4:6)/DU*TU];      % adim

% Define the n° of points for the Interpolation
opt.N = 1000;

% Interpolate the Ephemeris and Retrieve Target's Initial State
[X0t_MCI, COE0t, MEE0t, EarthPPsMCI, DSGPPsMCI, SunPPsMCI, MoonPPsECI, time, t0, tf, opt.N] = ...
    EphemerisHandler(deltaE, psiM, deltaM, opt.N, date0, datef);

% Combine the Target and Chaser States into TC ~ [Target; Chaser] State
TC0 = [MEE0t; RHO0_LVLH];

% Open log file
log = fopen('Output/log.txt', 'w+');

% Define the timespan for the propagation
tspan_ref = linspace(t0, tf, opt.N);
dt_ref = tspan_ref(2) - tspan_ref(1);


%% Propagate Reference Target Trajectory

% Perform the Propagation of the Target Trajectory
pbar = waitbar(0, 'Performing the Target Trajectory Propagation');
[~, MEEt_ref] = ode113(@(t, MEE) DynamicalModelMEE(t, MEE, EarthPPsMCI, SunPPsMCI, ...
    muE, muS, MoonPPsECI, deltaE, psiM, deltaM, t0, tf), tspan_ref, MEE0t, OptionsODE);
close(pbar);

% Conversion from MEE to COE
COEt_ref = MEE2COE(MEEt_ref);

% Conversion from COE to MCI
Xt_MCI_ref = COE2rvPCI(COEt_ref, muM);

% Interpolate the Angular Velocity of LVLH wrt MCI
[~, omegadotPPsLVLH] = TargetHandler(Xt_MCI_ref, COEt_ref, MEEt_ref, tspan_ref, ...
    EarthPPsMCI, SunPPsMCI, MoonPPsECI, deltaE, psiM, deltaM, muE, muS);


%% Chaser State Backwards Propagation from Final Conditions

% Define the Backwards Drift tspan
tspan_back = linspace(tf, t0 + (tf-t0)/2, opt.N);

% Define Backward Drift Propagation Settings
optODE_back = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'Events', @(t, Y) driftstop_back(t, Y));

% Define Final Conditions
MEEft = MEEt_ref(end, :)';
TCf_backdrift = [MEEft; RHOf_LVLH];

% Perform the Backpropagation of Target and Chaser Trajectories
[tspan_back, TC_backdrift] = ode113(@(t, TC) NaturalRelativeMotion(t, TC, EarthPPsMCI, SunPPsMCI, ...
    muE, muS, MoonPPsECI, deltaE, psiM, deltaM, t0, tf, omegadotPPsLVLH), tspan_back, TCf_backdrift, optODE_back);

% Check Successful Back Drift
if length(tspan_back) >= opt.N
    fprintf('Could not complete Backwards Natural Drift with the given Initial Conditions.')
end

% Retrieve the States and epoch for Final Drift
TC0_backdrift = TC_backdrift(end, :)';
t0_backdrift = tspan_back(end);             % !!! there's a 2e-4 difference


%% Generate the Reference Chaser Trajectory

% Set the Via Points and Interpolate the Reference Trajectory
[RHOdPPsLVLH, viapoints, t_viapoints] = ChebyschevReferenceTrajectory(TC0, TC0_backdrift, t0, t0_backdrift);

% Compute the Desired State Derivative
rhodot_rPPs = fnder(RHOdPPsLVLH(1), 1);
rhodot_tPPs = fnder(RHOdPPsLVLH(2), 1);
rhodot_hPPs = fnder(RHOdPPsLVLH(3), 1);
RHOdPPsLVLH = [RHOdPPsLVLH; rhodot_rPPs; rhodot_tPPs; rhodot_hPPs];

% Control Settings
prediction_interval = 120;                  % s - one prediction every 2 minutes
prediction_dt = prediction_interval/TU;     % adim - prediction interval adimensionalized

prediction_maxlength = fix(seconds(total_time)/prediction_interval);        % max length of the predictive controls
check_times = zeros(prediction_maxlength, 1);           % gets all the times at which one should check in continue to saturate by future prediction

for k = 1 : prediction_maxlength
    check_times(k) = t0 + (k-1) * prediction_dt;        % compute the check_times vector
end

prediction_delta = (30*60)/TU;      % predictive propagation of 30 min forward
N_inner = round(prediction_delta/dt_ref) - 1;       % n° of points in the inner propagation

% Parallel Computing Utils
fileNum = 1;
parallel_path = ['Data/temp/VelocitySim' num2str(fileNum) '.mat'];


%% Perform Chaser Rendezvous Manoeuvre and Natural Drift

% Define Initial Conditions
x70 = 1;                % initial mass ratio
TCC0 = [TC0; x70];      % initial TargetControlledChaser State

% Perform Rendezvous Propagation
[tspan_ctrl, TCC_ctrl] = ode_Ham(@(t, state) HybridPredictiveControl(t, state, EarthPPsMCI, SunPPsMCI, muM, ...
    muE, muS, time, MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, RHOdPPsLVLH, DU, TU, check_times, prediction_delta, N_inner, parallel_path),...
    [t0, t0_backdrift], TCC0, opt.N);

% Retrieve Final Mass Ratio Value
x7f = TCC_ctrl(end, 13);
bookmark = length(tspan_ctrl);


% Define Backward Drift Propagation Settings
optODE_fwd = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'Events', @(t, Y) driftstop_fwd(t, Y));

% Define the Forward Drift tspan
tspan_drift = linspace(tspan_ctrl(end), tf + abs(tf-tspan_ctrl(end)), opt.N);

% Define Initial Conditions
TC0_drift = TCC_ctrl(end, 1:12);

% Perform the Final Natural Drift of the Chaser
[tspan_drift, TC_drift] = ode113(@(t, TC) NaturalRelativeMotion(t, TC, EarthPPsMCI, SunPPsMCI, ...
    muE, muS, MoonPPsECI, deltaE, psiM, deltaM, t0, tf, omegadotPPsLVLH), tspan_drift, TC0_drift, optODE_fwd);

% Check Successful Forward Drift
if length(tspan_drift) >= opt.N         % !!! qui dovrei controllare lo stato o l'uscita dall'evento
    fprintf('Could not complete Forward Natural Drift with the given Initial Conditions.')
end


% Stack the States
TCC_drift = [TC_drift, x7f*ones(length(tspan_drift), 1)];
TCC = [TCC_ctrl; TCC_drift];
tspan = [tspan_ctrl; tspan_drift];


save('Data/temp/Raw Propagation.mat');


%% Post-Processing

% Retrieve States from Propagation
MEEt = TCC(:, 1:6);
RHO_LVLH = TCC(:, 7:12);
x7 = TCC(:, 13);

% Retrieve Target MCI State
COEt = MEE2COE(MEEt);
Xt_MCI = COE2rvPCI(COEt, muM);

M = length(tspan);

% Initialize Post-Processing Variables
RHO_MCI = zeros(M, 6);
RHOd_LVLH = zeros(M, 6);
Xc_MCI = zeros(M, 6);
dist = zeros(M, 1);
vel = zeros(M, 1);
u = zeros(3, M);
u_norms = zeros(M, 1);
f = zeros(3, M);
f_norms = zeros(M, 1);
kp_store = zeros(1, M);
kp_type = zeros(1, M);

load(parallel_path)     % load stopSaturationTime

g0 = 9.80665;
u_limit = 5e-5*g0*ones(M, 1);

% Initialize Terminal Trajectory Variables
terminal_tol = 50e-3/DU;    % tolerance for terminal trajectory flag
RHO_LVLH_T = [];
RHOd_LVLH_T = [];
terminal_flag = zeros(M, 1);
terminal_signchange = 0;

% Perform Post-Processing
pbar = waitbar(0, 'Performing the Final Post-Processing');
for i = 1 : size(RHO_LVLH, 1)
    
    % Compute MCI States
    RHO_MCI(i, :) = rhoLVLH2MCI(RHO_LVLH(i, :)', Xt_MCI(i, :)', tspan(i), EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM)';
    Xc_MCI(i, :) = Xt_MCI(i, :) + RHO_MCI(i, :);
    
    % Compute the Distance and Velocity norms
    dist(i) = norm(RHO_MCI(i, 1:3));
    vel(i) = norm(RHO_MCI(i, 4:6));

    % Compute the Desired State from Via Points interpolation
    RHOd_LVLH(i, :) = ppsval(RHOdPPsLVLH, tspan(i));

    % Control Quantities
    if i <= bookmark        % hybrid-predictive control law

        [~, ~, ~, ~, ~, ~, u(:, i), ~, ~, ~, f, f_norms(i), kp_store(i), kp_type(i)] = ...
            HybridPredictivePostProcessing(tspan_ctrl(i), TCC_ctrl(i, :), EarthPPsMCI, SunPPsMCI, muM, ...
            muE, muS, time, MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, RHOdPPsLVLH, DU, TU, check_times, stopSaturationTime);

        u_norms(i) = norm(u(:, i));

    elseif i > bookmark     % natural control law

        s = i - bookmark;
        [~, ~, ~, ~, ~, ~, f(:, i)] = NaturalRelativeMotion_PostProcessing(tspan_drift(s), TC_drift(s,:)', EarthPPsMCI, SunPPsMCI, muE, muS, ...
            MoonPPsECI, deltaE, psiM, deltaM, t0, tf, omegadotPPsLVLH);
            
        f_norms(i) = norm(f(:, i));

    end

    % Emergency Sphere Crossing Check
    if dist(i) <= terminal_tol
        RHO_LVLH_T = [RHO_LVLH_T; RHO_LVLH(i, :)];
        RHOd_LVLH_T = [RHOd_LVLH_T; RHOd_LVLH(i, :)];
        terminal_flag(i) = 1;
        if i > 1
            if terminal_flag(i-1) == 0
                terminal_signchange = terminal_signchange + 1;
            end
        end
    else
        terminal_flag(i) = 0;
        if i > 1
            if terminal_flag(i-1) == 1
                terminal_signchange = terminal_signchange + 1;
            end
        end
    end
    
    % Update Progress Bar
    waitbarMessage = sprintf('Post-Processing Progress: %.2f%%\n', i/M*100);
    waitbar(i/M, pbar, waitbarMessage);      % update the waitbar

end
close(pbar)

fclose(log);


save('Data/temp/Post-Processed Propagation.mat');


%% Visualize the Results

close all
clear
clc

load('Data/temp/Post-Processed Propagation.mat');


% Draw the Target, Chaser and Reference Chaser Trajectories in MCI
figure('name', 'Trajectory in MCI Space')
T = DrawTrajMCI3D(Xt_MCI(:, 1:3)*DU, '#d1d1d1', '-.');
C = DrawTrajMCI3D(Xc_MCI(:, 1:3)*DU);
legend([T, C], {'Target Trajectory', 'Chaser Trajectory'}, 'location', 'best');
view([140, 30]);
title('Target and Chaser MCI Trajectories')
if opt.saveplots
    saveas(gcf, strcat('Output/Trajectory MCI.jpg'))
end


% Visualize Chaser State in LVLH
figure('name', 'Chaser Trajectory in LVLH Space')
C_LVLH = DrawTrajLVLH3D(RHO_LVLH(:, 1:3)*DU);
Cd_LVLH = DrawTrajLVLH3D(RHOd_LVLH(:, 1:3)*DU, '#6efad2', '-.');
title('Chaser LVLH Trajectory')
if opt.saveplots
    saveas(gcf, strcat('Output/Trajectory LVLH.jpg'))
end


% Visualize the Terminal Chaser State in LVLH
figure('name', 'Terminal Chaser Trajectory in LVLH Space')
C_LVLH_T = DrawTrajLVLH3D(RHO_LVLH_T(:, 1:3)*DU);
Cd_LVLH_T = DrawTrajLVLH3D(RHOd_LVLH_T(:, 1:3)*DU, '#6efad2', '-.');
title('Terminal Chaser LVLH Trajectory')
if opt.saveplots
    saveas(gcf, strcat('Output/Trajectory Terminal LVLH.jpg'))
end


% Visualize LVLH State Components
figure('name', 'Chaser LVLH State Components')
subplot(2, 3, 1)
plot((tspan-t0)*TU*sec2hrs, RHOd_LVLH(:, 1)*DU, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
hold on
plot((tspan-t0)*TU*sec2hrs, RHO_LVLH(:, 1)*DU, 'color', '#4195e8', 'LineWidth', 1.5)
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\rho_r \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('Desired', 'Actual', 'fontsize', 10, 'location', 'best')

subplot(2, 3, 2)
plot((tspan-t0)*TU*sec2hrs, RHOd_LVLH(:, 2)*DU, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
hold on
plot((tspan-t0)*TU*sec2hrs, RHO_LVLH(:, 2)*DU, 'color', '#4195e8', 'LineWidth', 1.5)
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\rho_\theta \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('Desired', 'Actual', 'fontsize', 10, 'location', 'best')

subplot(2, 3, 3)
plot((tspan-t0)*TU*sec2hrs, RHOd_LVLH(:, 3)*DU, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
hold on
plot((tspan-t0)*TU*sec2hrs, RHO_LVLH(:, 3)*DU, 'color', '#4195e8', 'LineWidth', 1.5)
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\rho_h \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('Desired', 'Actual', 'fontsize', 10, 'location', 'best')

subplot(2, 3, 4)
plot((tspan-t0)*TU*sec2hrs, RHOd_LVLH(:, 4)*DU/TU, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
hold on
plot((tspan-t0)*TU*sec2hrs, RHO_LVLH(:, 4)*DU/TU, 'color', '#4195e8', 'LineWidth', 1.5)
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\dot{\rho}_r \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('Desired', 'Actual', 'fontsize', 10, 'location', 'best')

subplot(2, 3, 5)
plot((tspan-t0)*TU*sec2hrs, RHOd_LVLH(:, 5)*DU/TU, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
hold on
plot((tspan-t0)*TU*sec2hrs, RHO_LVLH(:, 5)*DU/TU, 'color', '#4195e8', 'LineWidth', 1.5)
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\dot{\rho}_\theta \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('Desired', 'Actual', 'fontsize', 10, 'location', 'best')

subplot(2, 3, 6)
plot((tspan-t0)*TU*sec2hrs, RHOd_LVLH(:, 6)*DU/TU, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
hold on
plot((tspan-t0)*TU*sec2hrs, RHO_LVLH(:, 6)*DU/TU, 'color', '#4195e8', 'LineWidth', 1.5)
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\dot{\rho}_h \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('Desired', 'Actual', 'fontsize', 10, 'location', 'best')
if opt.saveplots
    saveas(gcf, strcat('Output/State LVLH Components.jpg'))
end


% Draw the norms of distance and velocity
figure('Name', 'Distance and Velocity Norms')
subplot(1, 2, 1)
plot((tspan-t0)*TU*sec2hrs, dist*DU, 'color', '#4195e8', 'LineWidth', 1.5)
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$|\rho| \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
title('Relative Distance')
subplot(1, 2, 2)
plot((tspan-t0)*TU*sec2hrs, vel*DU/TU, 'color', '#4195e8', 'LineWidth', 1.5)
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$|v| \ [km/s]$', 'interpreter', 'latex', 'fontsize', 12)
title('Relative Velocity')
if opt.saveplots
    saveas(gcf, strcat('Output/Relative Distance and Velocity.jpg'))
end


% Visualize Control Norm
figure('name', 'Control Thrust')
p1 = plot((tspan-t0)*TU*sec2hrs, u_norms*1000*DU/TU^2, 'Color', '#4195e8', 'LineWidth', 1.5);
hold on
p2 = plot((tspan-t0)*TU*sec2hrs, u_limit, 'r--', 'LineWidth', 1.2);
p3 = plot((tspan-t0)*TU*sec2hrs, f_norms*1000*DU/TU^2, 'Color', '#93faad', 'LineWidth', 1.5);
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$[m/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
title('Control Norm')
legend([p1, p2, p3], '$|u|$', '$u_{max}$', '$f$','Location', 'best', 'Fontsize', 12, 'Interpreter','latex');
if opt.saveplots
    saveas(gcf, strcat('Output/Control Norm.jpg'))
end


% Visualize Control Components
figure('name', 'Control Thrust Components')
u1 = plot((tspan-t0)*TU*sec2hrs, u(1,:)*1000*DU/TU^2, 'LineWidth', 1.5);
hold on
u2 = plot((tspan-t0)*TU*sec2hrs, u(2,:)*1000*DU/TU^2, 'LineWidth', 1.5);
u3 = plot((tspan-t0)*TU*sec2hrs, u(3,:)*1000*DU/TU^2, 'LineWidth', 1.5);
ulim = plot((tspan-t0)*TU*sec2hrs, u_limit, 'r--', 'LineWidth', 1.2);
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$[m/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
title('Control Components')
legend([u1, u2, u3, ulim], '$u_r$', '$u_{\theta}$', '$u_h$', '$u_{max}$','Location', 'best', 'Fontsize', 12, 'Interpreter','latex');
if opt.saveplots
    saveas(gcf, strcat('Output/Control Components.jpg'))
end





if opt.create_animation
    figure('name', 'Rendezvous and Docking Animation', 'WindowState', 'maximized')
    DrawRendezvous(Xt_MCI(:, 1:3)*DU, Xc_MCI(:, 1:3)*DU, RHO_LVLH, bookmark, opt)
end


