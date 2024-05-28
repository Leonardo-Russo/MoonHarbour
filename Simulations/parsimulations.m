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


% Define the n° of simulations
MC = 100;

% Preallocate Data
data = struct('RHO_LVLH', [], 'M_ctrl_DA', [], 'DU', [], 'RHOd_LVLH', [], 'color', [], 'dist', [], 'vel', [], 'successful', [], 'deltaState', [], 'failure_times', [], 'misalignments', []);
data = repmat(data, MC, 1);

% Preallocate Table
table = zeros(MC, 10);

% Define colormap
cmap = lines(MC);

% Success Flags
successful_dist_tol = 1e-1;     % 10 cm
successful_vel_tol = 1e-2;      % 1 cm/s
success_color = '#39e617';
failure_color = '#d12828';

% % Parallel Computing
% pool = gcp('nocreate');
% if isempty(pool)
%     c = parcluster;                         % get the current cluster object
%     numDesiredWorkers = c.NumWorkers;       % get the maximum number of workers
%     pool = parpool(numDesiredWorkers);      % start the parallel pool with the maximum number of workers   
% end

% Define Simulation Options
sim_id = "combined_10s";
mkdir(strcat("Results/", sim_id));
sampling_time = 10;                     % seconds
verbose = true;
misalignment_type = "oscillating";
state_perturbation_flag = true;
engine_failure_flag = true;
include_actuation = false;

% parfor (mc = 1 : MC, pool.NumWorkers)
for mc = 1 : MC

    % try

        [RHO_LVLH, M_ctrl_DA, M_ctrl, DU, RHOd_LVLH, dist, vel, deltaState, failure_times, misalignments] = parfmain(sampling_time, include_actuation, verbose, misalignment_type, state_perturbation_flag, engine_failure_flag);

        is_safe = check_min_distance(dist, DU, M_ctrl_DA, M_ctrl, 9.8);

        data(mc).color = cmap(mc, :);
        data(mc).RHO_LVLH = RHO_LVLH;
        data(mc).M_ctrl_DA = M_ctrl_DA;
        data(mc).DU = DU;
        data(mc).RHOd_LVLH = RHOd_LVLH;
        data(mc).dist = dist;
        data(mc).vel = vel;
    
        if norm(deltaState(1:3)) <= successful_dist_tol && norm(deltaState(4:6)) <= successful_vel_tol
            data(mc).status = true;
        else
            data(mc).status = false;
        end
    
        if ~is_safe
            data(mc).status = data(mc).status - 0.5;
        end
    
        table(mc, :) = [mc, data(mc).status, deltaState, norm(deltaState(1:3)), norm(deltaState(4:6))];

    % catch
    % 
    %     fprintf('Simulation n° %2d was not successful.\n', mc);
    %     data(mc).status = -1;
    %     table(mc, :) = [mc, data(mc).status, zeros(1, 8)];
    % 
    % end


end

save(strcat("Results/", sim_id, "/", sim_id, ".mat"));
