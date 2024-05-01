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
colors = lines(MC);     % alternatively 'lines', 'jet', 'hsv', ...

parfor mc = 1 : MC

    workspace_path = strcat("Results/", sim_dir, "/", sim_dir, "-", string(mc), ".mat");

    fmain(workspace_path, false);

end