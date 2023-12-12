%% Modified Equinoctial Elements Orbital Propagator - Leonardo Russo

close all
clear
clc

addpath('Library/')
addpath('Data/')
addpath('Data/Planets/')

% % Set path to Python in the Specific Conda Enviroment
% pyversion('C:\Users\leo31\anaconda3\envs\MachineLearning310\python.exe')

%% Parameters and Settings Definition

% Define options for ode113()
OptionsODE = odeset('RelTol', 1e-7, 'AbsTol', 1e-6, 'MaxStep', Inf);

savechoice = 0;     % set as 1 to save a copy of the plots locally

% Define Global Variables
global DU TU Rm muM

% Define HyperParameters
DU = 10000;                                             % km
TU = sqrt(1 / 4902.7779 * DU^3);                        % s
Rm = 1738.1 / DU;                                       % DU
muM = 4902.7779 * TU^2 / DU^3;                          % DU^3/TU^2
muE = 398600.4418 * TU^2 / DU^3;                        % DU^3/TU^2
muS = 132712440018 * TU^2 / DU^3;                       % DU^3/TU^2
deltaE = deg2rad(23.45);                                % rad
psiM = deg2rad(-81.7 + 360/18.6 * (5 + 4.5/12));        % rad
deltaM = deg2rad(1.5);                                  % rad

Day = 86400;                                            % s

%% Ephemeris Interpolation

% Define the n° of points for the Interpolation
Npoints = 100000;

% Handle the Interpolation and Retrieve IC as MEE0
[MEE0, EarthPPsMCI, DSGPPsMCI, SunPPsMCI, MoonPPsECI, time, t0, tf] = ...
 EphemerisHandler(deltaE, psiM, deltaM, Npoints);

%% Propagation of the Trajectory

% Define the timespan for the propagation
tspan = linspace(t0, tf, Npoints);

% Create Progress Bar
global pbar
pbar = waitbar(0, 'Performing the Propagation');

% Perform the Integration
[~, MEE] = ode45(@(t, MEE) DynamicalModelMEE(t, MEE, EarthPPsMCI, SunPPsMCI, ...
                                             muE, muS, time, MoonPPsECI, deltaE, ...
                                             psiM, deltaM, t0, tf), tspan, MEE0, OptionsODE);
close(pbar);

%% Processing of the Results

% Conversion from MEE to COE
COE = MEE2COE(MEE);

% Conversion from COE to MCI
X_MCI = COE2rvPCI(COE, muM);

% Compute State in Synodic Reference Frame
rM2DSG_SYN = MCI2SYN(X_MCI, tspan, MoonPPsECI, deltaE, psiM, deltaM);

% Retrieve Reference DSG state in MCI - ppsval()
XDSG_MCI = zeros(length(tspan), 6);
for i = 1 : length(tspan)
    XDSG_MCI(i, :) = ppsval(DSGPPsMCI, tspan(i));
end

% Compute Perturbing Accelerations for Visualization
a3B_E = zeros(length(tspan), 3);
a3B_S = zeros(length(tspan), 3);
aG_M = zeros(length(tspan), 3);
normEs = zeros(length(tspan), 1);
normSs = zeros(length(tspan), 1);
normMs = zeros(length(tspan), 1);

for i = 1 : length(tspan)
    a3B_E(i, :) = ThirdBody(MEE(i, :), tspan(i), EarthPPsMCI, muE, tspan)'*DU/TU^2;                              % Third body
    a3B_S(i, :) = ThirdBody(MEE(i, :), tspan(i), SunPPsMCI, muS, tspan)'*DU/TU^2;                                % Fourth body
    aG_M(i, :) = MoonHarmPerts(MEE(i, :)', MoonPPsECI, tspan(i), tspan, muM, deltaE, psiM, deltaM)'*DU/TU^2;      % Moon Harmonics
    
    normEs(i) = norm(a3B_E(i, :));
    normSs(i) = norm(a3B_S(i, :));
    normMs(i) = norm(aG_M(i, :));
end

%% Visualization of the Results

close all
clc

% Load the Reference Trajectory into the current workspace
load('Reference.mat', 'Reference_State');   % recall CartState
DU_ref = 10000;
TU_ref = 14281.6665897583;


% Visualize Trajectory in MCI Reference Frame
figure('name', 'Trajectory in MCI Space')
P = DrawTrajMCI3D(X_MCI(:, 1:3)*DU);
R = DrawTrajMCI3D(Reference_State(:, 1:3)*DU_ref, '#3094cf', '--');
E = DrawTrajMCI3D(XDSG_MCI(:, 1:3)*DU, '#d1d1d1', ':');
legend([P, R, E], {'Propagated', 'Reference', 'Ephemeris'});
if savechoice
    saveas(gcf, strcat('Output/TrajMCI.jpg'))
end

final_distanceEPH = norm(X_MCI(end, 1:3)-XDSG_MCI(end, 1:3))*DU;
fprintf('Final Distance from Ephemeris is: %.2f km\n', final_distanceEPH);


% Visualize Trajectory in Synodic Reference Frame
figure('name', 'Trajectory in SYN Space')
S = DrawTrajMCI3D(rM2DSG_SYN(:, 1:3)*DU);
legend([S], {'Propagated'});
if savechoice
    saveas(gcf, strcat('Output/TrajSYN.jpg'))
end


% Visualize MCI State Errors wrt Ephemeris
figure('name', 'Carthesian State Error wrt Ephemeris')
subplot(2, 1, 1)
plot((tspan - t0) * TU / Day, (X_MCI(:, 1) - XDSG_MCI(:, 1)) * DU)
hold on
grid on
plot((tspan - t0) * TU / Day, (X_MCI(:, 2) - XDSG_MCI(:, 2)) * DU)
plot((tspan - t0) * TU / Day, (X_MCI(:, 3) - XDSG_MCI(:, 3)) * DU)
xlabel('t [Days]')
ylabel('[km]')
legend('$\delta_x$', '$\delta_y$', '$\delta_z$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')

subplot(2, 1, 2)
plot((tspan - t0) * TU / Day, (X_MCI(:, 4) - XDSG_MCI(:, 4)) * DU/TU)
hold on
grid on
plot((tspan - t0) * TU / Day, (X_MCI(:, 5) - XDSG_MCI(:, 5)) * DU/TU)
plot((tspan - t0) * TU / Day, (X_MCI(:, 6) - XDSG_MCI(:, 6)) * DU/TU)
xlabel('t [Days]')
ylabel('[km/s]')
legend('$\delta_{v_x}$', '$\delta_{v_y}$', '$\delta_{v_z}$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')
if savechoice
    saveas(gcf, strcat('Output/MCI_err_eph.jpg'))
end

% Visualize MCI State Errors wrt Reference
figure('name', 'Carthesian State Error wrt Reference')
subplot(2, 1, 1)
plot((tspan - t0) * TU / Day, X_MCI(:, 1)*DU - Reference_State(:, 1)*DU_ref)
hold on
grid on
plot((tspan - t0) * TU / Day, X_MCI(:, 2)*DU - Reference_State(:, 2)*DU_ref)
plot((tspan - t0) * TU / Day, X_MCI(:, 3)*DU - Reference_State(:, 3)*DU_ref)
xlabel('t [Days]')
ylabel('[km]')
legend('$\delta_x$', '$\delta_y$', '$\delta_z$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')

subplot(2, 1, 2)
plot((tspan - t0) * TU / Day, X_MCI(:, 4)*DU/TU - Reference_State(:, 4)*DU_ref/TU_ref)
hold on
grid on
plot((tspan - t0) * TU / Day, X_MCI(:, 5)*DU/TU - Reference_State(:, 5)*DU_ref/TU_ref)
plot((tspan - t0) * TU / Day, X_MCI(:, 6)*DU/TU - Reference_State(:, 6)*DU_ref/TU_ref)
xlabel('t [Days]')
ylabel('[km/s]')
legend('$\delta_{v_x}$', '$\delta_{v_y}$', '$\delta_{v_z}$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')
if savechoice
    saveas(gcf, strcat('Output/MCI_err_ref.jpg'))
end


% Visualize COE Errors wrt Ephemeris
figure('name', 'Evolution of COE wrt COE_REF')
COE_EPH = rvPCI2COE(XDSG_MCI, muM);
plotCOEvsREF(COE, COE_EPH, (tspan - t0) * TU / Day);
if savechoice
    saveas(gcf, strcat('Output/COEState_error.jpg'))
end


% Visualize MEE Errors wrt Errors
figure('name', 'Evolution of MEE wrt MEE_REF')
MEE_EPH = COE2MEE(COE_EPH);
plotMEEvsREF(MEE, MEE_EPH, (tspan - t0) * TU / Day);
if savechoice
    saveas(gcf, strcat('Output/MEEState_error.jpg'))
end


% Visualize the Perturbations
figure('name', 'Perturbation Accelerations')
subplot(2, 2, 1)
plot((tspan - t0) * TU / Day, a3B_E)
xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 11)
ylabel('$a_i \ [km/s^2]$', 'Interpreter','latex', 'FontSize', 11)
legend('$a_{pE_r}$', '$a_{pE_{\theta}}$', '$a_{pE_h}$', 'interpreter', 'latex', 'fontsize', 10, 'location', 'best')
title('Earth Perturbations')

subplot(2, 2, 2)
plot((tspan - t0) * TU / Day, a3B_S)
xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 11)
ylabel('$a_i \ [km/s^2]$', 'Interpreter','latex', 'FontSize', 11)
legend('$a_{pS_r}$', '$a_{pS_{\theta}}$', '$a_{pS_h}$', 'interpreter', 'latex', 'fontsize', 10, 'location', 'best')
title('Sun Perturbations')

subplot(2, 2, 3)
plot((tspan - t0) * TU / Day, aG_M)
xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 11)
ylabel('$a_i \ [km/s^2]$', 'Interpreter','latex', 'FontSize', 11)
legend('$a_{pM_r}$', '$a_{pM_{\theta}}$', '$a_{pM_h}$', 'interpreter', 'latex', 'fontsize', 10, 'location', 'best')
title('Moon Perturbations')

subplot(2, 2, 4)
plot((tspan - t0) * TU / Day, normEs)
hold on
plot((tspan - t0) * TU / Day, normSs)
plot((tspan - t0) * TU / Day, normMs)
xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 11)
ylabel('$|a_i| \ [km/s^2]$', 'Interpreter','latex', 'FontSize', 11)
title('Perturbation Norms')
legend('$|a_{pE}|$', '$|a_{pS}|$', '$|a_{pM}|$', 'interpreter', 'latex', 'fontsize', 10, 'location', 'best')
if savechoice
    saveas(gcf, strcat('Output/Perturbations.jpg'))
end


