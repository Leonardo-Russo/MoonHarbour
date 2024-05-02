%% MoonHarbour Simulations - Leonardo Russo

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

% Define the n° of simulations
MC = 10;

% Define Simulation Options
sim_dir = "misalignment";
mkdir(strcat("Results/", sim_dir));

parfor mc = 1 : MC

    workspace_path = strcat("Results/", sim_dir, "/", sim_dir, "-", string(mc), ".mat");

    fmain(workspace_path, false, "const", false, false);

    disp(mc)

end


% Define Simulation Options
sim_dir = "state_perturbation";
mkdir(strcat("Results/", sim_dir));

parfor mc = 1 : MC

    workspace_path = strcat("Results/", sim_dir, "/", sim_dir, "-", string(mc), ".mat");

    fmain(workspace_path, false, "null", true, false);

    disp(mc)

end



% Define Simulation Options
sim_dir = "engine_failure";
mkdir(strcat("Results/", sim_dir));

parfor mc = 1 : MC

    workspace_path = strcat("Results/", sim_dir, "/", sim_dir, "-", string(mc), ".mat");

    fmain(workspace_path, false, "null", false, true);

    disp(mc)

end