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
opt.additional_plots = false;
opt.showgui = false;
opt.N = 1000;                   % n° of points for the Interpolation
opt.RelTolODE = 1e-7;           % options for ode()
opt.AbsTolODE = 1e-6;


% Define the n° of simulations
MC = 100;

% Define Simulation Options
sim_dir = "state_perturbation";
mkdir(strcat("Results/", sim_dir));

for mc = 1 : MC

    workspace_path = strcat("Results/", sim_dir, "/", sim_dir, "-", string(mc), ".mat");

    fmain(workspace_path, false, "null", true, false);

    disp(mc);

end


% Define Simulation Options
sim_dir = "engine_failure";
mkdir(strcat("Results/", sim_dir));

for mc = 1 : MC

    workspace_path = strcat("Results/", sim_dir, "/", sim_dir, "-", string(mc), ".mat");

    fmain(workspace_path, false, "null", false, true);

    disp(mc);

end
