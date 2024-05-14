%% MoonHarbour Simulations - Leonardo Russo

close all
clear
clc
                      
addpath('../')
addpath('../Library/')
addpath('../Data/')
addpath('../Data/Planets/')
addpath('../Data/Materials/')
addpath('../Data/Ephemeris/')
addpath('../Data/Utils/')


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


% Define the n° of simulations
MC = 10;

% numDesiredWorkers = 4; % Example: set to desired number of workers
% pool = parpool(numDesiredWorkers); % Starts a pool with the specified number of workers

% Define Simulation Options
sim_dir = "misalignment_15s";
mkdir(strcat("Results/", sim_dir));
sampling_time = 15;                     % seconds
verbose = true;
misalignment_type = "oscillating";
state_perturbation_flag = false;
engine_failure_flag = false;

for i = 2 : MC

    workspace_path = strcat("Results/", sim_dir, "/", sim_dir, "-", string(i), ".mat");

    fmain(sampling_time, workspace_path, verbose, misalignment_type, state_perturbation_flag, engine_failure_flag);

    disp(i);

end
