%% Relative Motion Propagator - Leonardo Russo

close all
clear
clc

addpath('Library/')
addpath('Data/')
addpath('Data/Planets/')
addpath('Data/Materials/')
addpath('Dario/')
addpath('temp/')

skip = 1;

if ~skip

    % Introduce Options Structure
    opt = struct('name', "Progator Options");
    opt.saveplots = false;
    
    % Define options for ode113()
    opt.RelTolODE = 1e-7;
    opt.AbsTolODE = 1e-6;
    OptionsODE = odeset('RelTol', opt.RelTolODE, 'AbsTol', opt.AbsTolODE);
    
    
    %% Hyperparameters and Settings
    
    % Define Global Variables
    global DU TU Rm muM pbar log
    
    % Define Constants
    DU = 10000;                                             % km
    TU = sqrt(1 / 4902.7779 * DU^3);                        % s
    Rm = 1738.1 / DU;                                       % DU
    muM = 4902.7779 * TU^2 / DU^3;                          % DU^3/TU^2
    muE = 398600.4418 * TU^2 / DU^3;                        % DU^3/TU^2
    muS = 132712440018 * TU^2 / DU^3;                       % DU^3/TU^2
    deltaE = deg2rad(23.45);                                % rad
    psiM = deg2rad(-81.7 + 360/18.6 * (5 + 4.5/12));        % rad
    deltaM = deg2rad(1.5);                                  % rad
    Dsol = 86400;                                           % s
    
    % Define Time Domain
    date0 = datetime('2025-05-23 9:56:10'); 
    datef = datetime('2025-05-23 21:56:10');    % @ periselenium
    total_time = datef - date0;
    
    % Define the Chaser Initial Conditions
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
    
    % Combine the Target and Chaser States into Y State
    Y0 = [MEE0t; RHO0_LVLH];
    
    % Open log file
    log = fopen('Data/log.txt', 'w+');
    
    % Define the timespan for the propagation
    tspan = linspace(t0, tf, opt.N);
    
    
    %% Propagate Target Trajectory using MEE
    
    % Perform the Propagation of the Target Trajectory
    pbar = waitbar(0, 'Performing the Target Trajectory Propagation');
    [~, MEEt] = ode113(@(t, MEE) DynamicalModelMEE(t, MEE, EarthPPsMCI, SunPPsMCI, ...
    muE, muS, MoonPPsECI, deltaE, psiM, deltaM, t0, tf), tspan, MEE0t, OptionsODE);
    close(pbar);
    
    % Conversion from MEE to COE
    COEt = MEE2COE(MEEt);
    
    % Conversion from COE to MCI
    Xt_MCI = COE2rvPCI(COEt, muM);
    
    % Interpolate the Angular Velocity of LVLH wrt MCI
    [omegaPPsLVLH, omegadotPPsLVLH] = TargetHandler(Xt_MCI, COEt, MEEt, tspan, ...
    EarthPPsMCI, SunPPsMCI, MoonPPsECI, deltaE, psiM, deltaM, muE, muS);
    
    % Visualize Target Trajectory during R&D
    figure('name', 'Target MCI Trajectory during Rendezvous and Docking')
    T = DrawTrajMCI3D(Xt_MCI(:, 1:3)*DU);
    
    
    %% Backwards Propagation from Final Conditions
    
    % Define the Backwards tspan
    tspan_back = linspace(tf, t0 + (tf-t0)/2, opt.N);
    
    % Define Backward Propagation Settings
    optODE_back = odeset('RelTol', opt.RelTolODE, 'AbsTol', opt.AbsTolODE, 'Events', @(t, Y) backstop(t, Y));
    
    % Define Final Conditions
    MEEft = MEEt(end, :)';
    Yf = [MEEft; RHOf_LVLH];
    
    % Perform the Backpropagation of Target and Chaser Trajectories
    pbar = waitbar(0, 'Performing the Chaser Trajectory Propagation');
    [tspan_back, Y] = ode113(@(t, Y) NaturalRelativeMotion(t, Y, EarthPPsMCI, SunPPsMCI, ...
    muE, muS, MoonPPsECI, deltaE, psiM, deltaM, t0, tf, omegadotPPsLVLH), tspan_back, Yf, optODE_back);
    close(pbar);
    
    % Retrieve the States and epoch for Final Drift
    Ydrift = Y(end, :)';
    tdrift = tspan_back(end);
    
    % SHOW THE EVOLUTION IN THE DRIFT -> CREATE A DOCKING ANIMATION BASICALLY
    save("First Part.mat");

else

    load("First Part.mat");

end


%% Generate the Reference Trajectory


% Set the Via Points and Interpolate the Reference Trajectory
[ppXd, ViaPoints, t_via] = ChebyschevReferenceTrajectory(Y0, Ydrift, t0, tdrift);

% Define Forward Propagation Settings
optODE_fwd = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

% Control Settings
prediction_step = 120;              % prediction interval in seconds -> one prediction every 2 minutes
Delta_t = prediction_step/TU;
len0 = fix(seconds(total_time)/prediction_step);      % length of the predictive controls
checkTimes = zeros(len0, 1);        % gets all the times at which one should check in continue to saturate by future prediction

for p = 1 : len0
    checkTimes(p) = t0 + (p-1) * Delta_t;
end

Delta_integration = (30*60)/TU; % predictive propagation of 30 min

N_inner_integration = round(Delta_integration/(tspan(2)-tspan(1))) - 1;

fileNum = 1;
P = ['Dario\VelocitySim'...
        num2str(fileNum) '.mat'];
filedir = 'Dario\';
filename = ['VelocitySim' num2str(fileNum) '.mat'];

%% g

x70 = 1;            % initial mass ratio
S0 = [Y0; x70];

[tspan, ControlledRelativeState] = ode_Ham(@(t, state) HybridPredictiveControl(t, state, EarthPPsMCI, SunPPsMCI, muM, ...
     muE, muS, time, MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, ppXd, DU, TU, checkTimes, Delta_integration, N_inner_integration, P),...
     [t0, tdrift], S0, opt.N);

tspan_control = tspan;


% Final natural propagation to see the real drift obtained
initialNaturalState = ControlledRelativeState(end,1:12);
options_nat2 = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'Events', @(t,state)stopper2(t,state,DU));
tspan_nat2 = linspace(tspan_control(end),tf, Npoints/10);
[t_drift, NaturalRelativeState2] = ode113(@(t, state) NaturalRelativeMotion(t, state, ppEarthMCI, ppSunMCI, muM, ...
    muE, muS, time, ppMoonECI, deltaE, psiM, deltaM, ppOmegaVect), tspan_nat2,...
    initialNaturalState, options_nat2);




