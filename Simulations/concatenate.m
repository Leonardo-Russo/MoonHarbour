%% Concatenate MC Simulations - Leonardo Russo

close all
clear
clc

addpath('../')
addpath('../Library/')
addpath('../Data/')
addpath('../Data/Planets/')
addpath('../Data/Materials/')
addpath('../Data/Ephemeris/')

root_dir = "Results/Partial Results";       % root results folder
sim_id1 = "berthing_60s_new_part1";
sim_id2 = "berthing_60s_new_part2";        % specific results identifier
sim_id_tot = "test";
idx1 = 51;
idx2 = 51;

% Preallocate Table
table_tot = zeros(idx1+idx2, 12);

% Define colormap
cmap_tot = lines(idx1+idx2);

% Load First Simulation
load(strcat(root_dir, "/", sim_id1, "/", sim_id1, ".mat"));
data_1 = data(1:idx1);
table_tot(1:idx1, :) = table(1:idx1, :);

clear data

% Load First Simulation
load(strcat(root_dir, "/", sim_id2, "/", sim_id2, ".mat"));
data_2 = data(1:idx2);
table_tot(idx1+1:end, :) = table(1:idx2, :);
table_tot(idx1+1:end, 1) = table_tot(idx1+1:end, 1) + idx1;

sim_id = sim_id_tot;
MC = idx1 + idx2;

% Concatenate the Results
table = table_tot;
data = [data_1; data_2];

% If you need to substitute stuff...
data(27) = data(101);
table(27, 2:end) = table(101, 2:end);
data(51+30) = data(102);
table(51+30, 2:end) = table(102, 2:end);
data = data(1:100);
table = table(1:100, :);
MC = 100;

for i = 1 : MC
    data(i).color = cmap_tot(i, :);
end

mkdir(strcat(root_dir, "/", sim_id_tot));
save(strcat(root_dir, "/", sim_id_tot, "/", sim_id_tot, ".mat"), 'sim_id', 'sampling_time', 'final_velocity', 'misalignment_type', 'engine_failure_flag', 'state_perturbation_flag', 'include_actuation', 'data', 'table', 'MC');






