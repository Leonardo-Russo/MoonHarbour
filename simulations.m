%% MoonHarbour Simulations - Leonardo Russo

close all
clear
clc
                                                                                
addpath('Library/')
addpath('Library/Dario/')
addpath('Data/')
addpath('Data/Planets/')
addpath('Data/Materials/')
addpath('Data/temp/')

% Define the n° of simulations
MC = 60;

% Generate Distributions
[beta_dist, gamma_dist] = define_distributions(MC);
colors = lines(MC);     % alternatively 'lines', 'jet', 'hsv', ...

% Initialize Stored Data
data = cell(MC, 1);

global DU

parfor mc = 1 : MC

    beta = beta_dist(mc);
    gamma = gamma_dist(mc);

    [RHO_LVLH, M_ctrl_DA] = fullsim(beta, gamma);
    
    data{mc} = {RHO_LVLH, M_ctrl_DA};

end

save('montecarlo.mat');


%% Visualize the Results

close all
clear
load('montecarlo.mat');
clc

saveplots = 1;
colors = hsv(MC);     % alternatively 'lines', 'jet', 'hsv', ...

% Visualize Chaser State in LVLH
figure('name', 'Chaser Trajectory in LVLH Space')
title('Chaser LVLH Trajectories')
for mc = 1 : MC
    RHO_LVLH = data{mc}{1};
    DrawTrajLVLH3D(RHO_LVLH(:, 1:3)*DU, colors(mc, :));    
end
if saveplots
    saveas(gcf, strcat('Output/Montecarlo - Trajectory LVLH.jpg'))
end


% Visualize the Terminal Chaser State in LVLH
figure('name', 'Terminal Chaser Trajectory in LVLH Space')
title('Terminal Chaser LVLH Trajectory')
for mc = 1 : MC
    RHO_LVLH = data{mc}{1};
    M_ctrl_DA = data{mc}{2};
    DrawTrajLVLH3D(RHO_LVLH(M_ctrl_DA:end, 1:3)*DU, colors(mc, :));
end
if saveplots
    saveas(gcf, strcat('Output/Montecarlo - Trajectory Terminal LVLH.jpg'))
end

