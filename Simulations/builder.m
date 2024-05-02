%% Builder for MC Simulations - Leonardo Russo

close all
clear
clc

addpath('../')
addpath('../Library/')
addpath('../Data/')
addpath('../Data/Planets/')
addpath('../Data/Materials/')
addpath('../Data/temp/')
addpath('../Simulations/')

root_dir = "Results";       % root results folder
sim_id = "engine_failure";    % specific results identifier


%% Extract the Results

% Get all files for the simulation
sim_dir = strcat(root_dir, "/", sim_id, "/");
files = dir(fullfile(sim_dir, '*.mat'));
files = files(~strcmp({files.name}, 'data.mat'));
n_sims = length(files);

% Preallocate Data
data = struct('RHO_LVLH', [], 'M_ctrl_DA', [], 'DU', [], 'RHOd_LVLH', [], 'color', [], 'dist', [], 'vel', [], 'successful', [], 'deltaState', [], 'failure_times', [], 'misalignments', []);
data = repmat(data, n_sims, 1);

% Preallocate Table
table = zeros(n_sims, 8);

% Define colormap
cmap = lines(n_sims);

successful_dist_tol = 1e-1;     % 10 cm
successful_vel_tol = 1e-2;      % 1 cm/s
success_color = '#39e617';
failure_color = '#d12828';


for k = 1 : n_sims

    temp = load(fullfile(files(k).folder, files(k).name), 'RHO_LVLH', 'M_ctrl_DA', 'DU', 'RHOd_LVLH', 'dist', 'vel', 'deltaState', 'failure_times', 'misalignments');
    
    data(k).color = cmap(k, :);
    data(k).RHO_LVLH = temp.RHO_LVLH;
    data(k).M_ctrl_DA = temp.M_ctrl_DA;
    data(k).DU = temp.DU;
    data(k).RHOd_LVLH = temp.RHOd_LVLH;
    data(k).dist = temp.dist;
    data(k).vel = temp.vel;

    if norm(temp.deltaState(1:3)) <= successful_dist_tol && norm(temp.deltaState(4:6)) <= successful_vel_tol
        data(k).successful = true;
    else
        data(k).successful = false;
    end

    table(k, :) = [k, data(k).successful, temp.deltaState];
    
end

save(strcat(sim_dir, "data.mat"));


%% Visualize the Simulations

% Show the table
disp(array2table(table, 'VariableNames', {'id', 'status', 'dr (m)', 'dtheta (m)', 'dh (m)', 'dv_r (m/s)', 'dv_theta (m/s)', 'dv_h (m/s)'}));

figure('name', 'Terminal Chaser Trajectory in LVLH Space')
title('Terminal Chaser LVLH Trajectory')
for k = 1 : n_sims
    if data(k).successful
        color = success_color;
    else
        color = failure_color;
    end
    DrawTrajLVLH3D(data(k).RHO_LVLH(data(k).M_ctrl_DA:end, 1:3) * data(k).DU, color);
end
view(-65, 15)
