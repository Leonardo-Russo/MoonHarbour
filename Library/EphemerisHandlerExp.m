function [X0_DSG, COE0, MEE0, EarthPPsMCI, DSGPPsMCI, SunPPsMCI, MoonPPsECI, time, t0,...
          tf, Npoints] = EphemerisHandlerExp(deltaE, psiM, deltaM, Npoints, date0, datef)
% Description: this function creates the necessary structures to continue
% witht the trajectory propagation. In order to reduce computational
% effort, it creates subvectors of Npoints from the original states and
% then continues the process with this subspace.
% 
% Inputs:
% deltaE = Earth's Ecliptic Obliquity
% psiM = Moon's Precession Angle
% deltaM = Moon's Equator Obliquity
% muE = Earth's gravitational parameter in canonical units
% muS = Sun's gravitational parameter in canonical units
% Npoints = n° of points for the subspace reduction
% 
% Outputs:
% MEE0 = initial value of the Modified Equinoctial Elements
% EarthPPsMCI = pp struct with the interpolation of the Earth State in MCI
% DSGPPsMCI = pp struct with the interpolation of the DSG State in MCI
% SunPPsMCI = pp struct with the interpolation of the Sun State in MCI
% MoonPPsECI = pp struct with the interpolation of the Moon State in ECI
% time = reduced time span
% t0 = initial time
% tf = final time

% Load Data from Ephemeris - states are provided in ECI
load('ephemeris_3M_17_27_05_2025.mat', 'stateDSG', 'stateMoon', 'stateSun', 'time');

% Recall Global Variables
global DU TU muM Rm

Dsol = 86400;       % s

% Look for Initial and Final times
t0_sharp = juliandate(date0)*Dsol;
tf_sharp = juliandate(datef)*Dsol;

t0 = interp1(time, time, t0_sharp, 'nearest');
t0_idx = find(time == t0);

tf = interp1(time, time, tf_sharp, 'nearest');
tf_idx = find(time == tf); 

tf = tf/TU;     % save the final time

tf_idx = tf_idx + (tf_idx - t0_idx) / 2;        % expand the final time index

% Retrieve Data from Ephemeris
stateDSG_ECI = stateDSG(t0_idx : tf_idx, :);
stateMoon_ECI = stateMoon(t0_idx : tf_idx, :);
stateSun_ECI = stateSun(t0_idx : tf_idx, :);
stateEarth_ECI = zeros(length(t0_idx : tf_idx), 6);
time = time(t0_idx : tf_idx);

% Set upper limit to Npoints
if Npoints > tf_idx - t0_idx
    Npoints = tf_idx - t0_idx;
    fprintf('Note! The n° of points for the interpolation exceeded \nthe n° of total points, therefore it has been set as:\nNpoints = %.0f\n', Npoints);
end

% Create a Data Subset 
indices = round(linspace(1, tf_idx - t0_idx, Npoints)');

time = time(indices);
stateDSG_ECI = stateDSG_ECI(indices, :);
stateMoon_ECI = stateMoon_ECI(indices, :);
stateSun_ECI = stateSun_ECI(indices, :);
stateEarth_ECI = stateEarth_ECI(indices, :);

% Conversion into Canonical Units
stateDSG_ECI(:, 1:3) = stateDSG_ECI(:, 1:3) / DU;
stateDSG_ECI(:, 4:6) = stateDSG_ECI(:, 4:6) / DU * TU;
stateMoon_ECI(:, 1:3) = stateMoon_ECI(:, 1:3) / DU;
stateMoon_ECI(:, 4:6) = stateMoon_ECI(:, 4:6) / DU * TU;
stateSun_ECI(:, 1:3) = stateSun_ECI(:, 1:3) / DU;
stateSun_ECI(:, 4:6) = stateSun_ECI(:, 4:6) / DU * TU;
stateEarth_ECI(:, 1:3) = stateEarth_ECI(:, 1:3) / DU;
stateEarth_ECI(:, 4:6) = stateEarth_ECI(:, 4:6) / DU * TU;

time = time / TU;
t0 = time(1);
% tf = time(end);

% Initialize Local Variables
stateEarth_MCI = zeros(Npoints, 6);
stateDSG_MCI = zeros(Npoints, 6);
stateSun_MCI = zeros(Npoints, 6);


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
        X0_DSG = [rDSG_MCI; vDSG_MCI];       % DSG State in MCI
        COE0 = rvPCI2COE(X0_DSG', muM);
        MEE0 = COE2MEE(COE0)';
    end
    

    % Assign values to the Matrices
    stateEarth_MCI(i, :) = [rE_MCI', vE_MCI'];
    stateSun_MCI(i, :) = [rS_MCI', vS_MCI'];
    stateDSG_MCI(i, :) = [rDSG_MCI', vDSG_MCI'];

end


% Perform the Interpolation
EarthPPsMCI = get_statePP(time, stateEarth_MCI);

SunPPsMCI = get_statePP(time, stateSun_MCI);

DSGPPsMCI = get_statePP(time, stateDSG_MCI);

MoonPPsECI = get_statePP(time, stateMoon_ECI);

end