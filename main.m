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

% Define the Chaser Initial Conditions
% Cerca di tenere questa condizione: norm(rhodot_LVLH) = 1 m/s -> il tuning
% dei parametri è stato fatto su questa ipotesi.
RHO0_LVLH = [1.5, 0, 0, 0, -1e-3, -1e-3]';                  % km, km/s
RHO0_LVLH = [RHO0_LVLH(1:3)/DU; RHO0_LVLH(4:6)/DU*TU];      % adim

% Define Desired Conditions for Docking
RHOf_LVLH = [5e-3, 0, 0, -1e-5, 0, 0]';                     % km, km/s
RHOf_LVLH = [RHOf_LVLH(1:3)/DU; RHOf_LVLH(4:6)/DU*TU];      % adim

% Define the n° of points for the Interpolation
opt.N = 1000;

% Interpolate the Ephemeris and Retrieve Target's Initial State
[X0t_MCI, COE0t, MEE0t, EarthPPsMCI, DSGPPsMCI, SunPPsMCI, MoonPPsECI, time, t0, tf, opt.N] = ...
    EphemerisHandler(deltaE, psiM, deltaM, opt.N, date0, datef);

% Combine the Target and Chaser States into TC ~ [Target; Chaser] State
TC0 = [MEE0t; RHO0_LVLH];

% Open log file
log = fopen('Data/log.txt', 'w+');

% Define the timespan for the propagation
tspan_ref = linspace(t0, tf, opt.N);
dt_ref = tspan_ref(2) - tspan_ref(1);


%% Propagate Reference Target Trajectory using MEE

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


%% Backwards Propagation from Final Conditions

% Define the Backwards Drift tspan
tspan_back = linspace(tf, t0 + (tf-t0)/2, opt.N);

% Define Backward Drift Propagation Settings
optODE_back = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'Events', @(t, Y) driftstop_back(t, Y));

% Define Final Conditions
MEEft = MEEt_ref(end, :)';
TCf_backdrift = [MEEft; RHOf_LVLH];

% Perform the Backpropagation of Target and Chaser Trajectories
pbar = waitbar(0, 'Performing the Chaser Trajectory Propagation');
[tspan_back, TC_backdrift] = ode113(@(t, TC) NaturalRelativeMotion(t, TC, EarthPPsMCI, SunPPsMCI, ...
    muE, muS, MoonPPsECI, deltaE, psiM, deltaM, t0, tf, omegadotPPsLVLH), tspan_back, TCf_backdrift, optODE_back);
close(pbar);

% Check Successful Back Drift
if length(tspan_back) >= opt.N
    error('Could not complete Backwards Natural Drift with the given Initial Conditions.')
end

% Retrieve the States and epoch for Final Drift
TC0_backdrift = TC_backdrift(end, :)';
t0_backdrift = tspan_back(end);             % !!! there's a 2e-4 difference


%% Generate the Reference Trajectory

% Set the Via Points and Interpolate the Reference Trajectory
[ppXd, ViaPoints, t_via] = ChebyschevReferenceTrajectory(TC0, TC0_backdrift, t0, t0_backdrift);

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


%% Perform Rendezvous Manoeuvre

% Define Initial Conditions
x70 = 1;                % initial mass ratio
TCC0 = [TC0; x70];      % initial TargetControlledChaser State

% Perform Rendezvous Propagation
[tspan_ctrl, TCC] = ode_Ham(@(t, state) HybridPredictiveControl(t, state, EarthPPsMCI, SunPPsMCI, muM, ...
    muE, muS, time, MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, ppXd, DU, TU, check_times, prediction_delta, N_inner, parallel_path),...
    [t0, t0_backdrift], TCC0, opt.N);

% Retrieve Final Mass Ratio Value
x7f = TCC(end, 13);
bookmark = length(tspan_ctrl);


% Define Backward Drift Propagation Settings
optODE_fwd = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'Events', @(t, Y) driftstop_fwd(t, Y));

% Define the Forward Drift tspan
tspan_drift = linspace(tspan_ctrl(end), tf + (tf-tspan_ctrl(end))/2, opt.N);

% Define Initial Conditions
TC0_drift = TCC(end, 1:12);

% Perform the Final Natural Drift of the Chaser
pbar = waitbar(0, 'Performing the Chaser Trajectory Propagation');
[tspan_drift, TC_drift] = ode113(@(t, TC) NaturalRelativeMotion(t, TC, EarthPPsMCI, SunPPsMCI, ...
    muE, muS, MoonPPsECI, deltaE, psiM, deltaM, t0, tf, omegadotPPsLVLH), tspan_drift, TC0_drift, optODE_fwd);
close(pbar);

% Check Successful Forward Drift
if length(tspan_drift) >= opt.N         % !!! qui dovrei controllare lo stato o l'uscita dall'evento
    error('Could not complete Forward Natural Drift with the given Initial Conditions.')
end

save('Data/temp/Propagation Completed.mat');


%% Post-Processing

% Stack the States
TCC_drift = [TC_drift, x7f*ones(length(tspan_drift), 1)];
TCC = [TCC; TCC_drift];
tspan = [tspan_ctrl; tspan_drift];

% Retrieve States from Propagation
MEEt = TCC(:, 1:6);
RHO_LVLH = TCC(:, 7:12);
x7 = TCC(:, 13);

% Retrieve Target MCI State
COEt = MEE2COE(MEEt);
Xt_MCI = COE2rvPCI(COEt, muM);

% Initialize Post-Processing Variables
RHO_MCI = zeros(length(tspan), 6);
Xc_MCI = zeros(length(tspan), 6);
dist = zeros(length(tspan), 1);

terminal_tol = 50e-3/DU;    % tolerance for terminal trajectory flag
RHOT_LVLH = [];
terminal_flag = zeros(length(tspan), 1);
terminal_signchange = 0;

% Perform Post-Processing
pbar = waitbar(0, 'Performing the Final Post-Processing');
for i = 1 : size(RHO_LVLH, 1)

    RHO_MCI(i, :) = rhoLVLH2MCI(RHO_LVLH(i, :)', Xt_MCI(i, :)', tspan(i), EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM)';
    Xc_MCI(i, :) = Xt_MCI(i, :) + RHO_MCI(i, :);

    dist(i) = norm(RHO_MCI(i, 1:3));

    if dist(i) <= terminal_tol
        RHOT_LVLH = [RHOT_LVLH; RHO_LVLH(i, :)];
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

    waitbarMessage = sprintf('Final Post-Processing Progress: %.2f%%\n', i/length(tspan)*100);
    waitbar(i/length(tspan), pbar, waitbarMessage);      % update the waitbar

end
close(pbar)

save('Data/temp/Post-Processing.mat');


%% Visualize the Results

close all
clear
clc

load('Data/temp/Post-Processing.mat');

if opt.create_animation
    figure('name', 'Rendezvous and Docking Animation', 'WindowState', 'maximized')
    DrawRendezvous(Xt_MCI(:, 1:3)*DU, Xc_MCI(:, 1:3)*DU, RHO_LVLH, bookmark, opt)
end


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


% Plot the evolution of the Chaser State in LVLH
figure('name', 'Chaser Trajectory in LVLH Space')
C_LVLH = DrawTrajLVLH3D(RHO_LVLH(:, 1:3)*DU);
title('Chaser LVLH Trajectory')
if opt.saveplots
    saveas(gcf, strcat('Output/Trajectory LVLH.jpg'))
end


% Plot the evolution of the Terminal Chaser State in LVLH
figure('name', 'Terminal Chaser Trajectory in LVLH Space')
C_LVLH = DrawTrajLVLH3D(RHOT_LVLH(:, 1:3)*DU);
title('Terminal Chaser LVLH Trajectory')
if opt.saveplots
    saveas(gcf, strcat('Output/Trajectory Terminal LVLH.jpg'))
end


% Visualize LVLH State
figure('name', 'Chaser LVLH State')
subplot(2, 1, 1)
plot((tspan - t0)*TU*sec2hrs, RHO_LVLH(:, 1)*DU)
hold on
grid on
plot((tspan - t0)*TU*sec2hrs, RHO_LVLH(:, 2)*DU)
plot((tspan - t0)*TU*sec2hrs, RHO_LVLH(:, 3)*DU)
title('Chaser LVLH Position')
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$[km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('$r$', '$\theta$', '$h$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')

subplot(2, 1, 2)
plot((tspan - t0)*TU*sec2hrs, RHO_LVLH(:, 4)*DU/TU)
hold on
grid on
plot((tspan - t0)*TU*sec2hrs, RHO_LVLH(:, 5)*DU/TU)
plot((tspan - t0)*TU*sec2hrs, RHO_LVLH(:, 6)*DU/TU)
title('Chaser LVLH Velocity')
xlabel('$t \ [days]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$[km/s]$', 'interpreter', 'latex', 'fontsize', 12)
legend('$\dot{r}$', '$\dot{\theta}$', '$\dot{h}$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')
if opt.saveplots
    saveas(gcf, strcat('Output/Chaser LVLH State.jpg'))
end





