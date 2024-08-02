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
MC = 32;

% Preallocate Data
data = struct('status', [], 'RHO_LVLH', [], 'M_ctrl_DA', [], 'DU', [], 'RHOd_LVLH', [], 'color', [], 'dist', [], 'vel', [], 'successful', [], 'deltaState', [], 'failure_times', [], 'misalignments', [], ...
              'renderdata', [], 'TCC', [], 'Xt_MCI', [], 'RHO_MCI', [], 'u', [], 'u_norms', [], 'f_norms', [], 'kp_store', [], 'qe0', [], 'qe', [], 'Tc', [], 'Ta', [], 'betas', [], 'gammas', [], 'acc', []);
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

% Parallel Computing
pool = gcp('nocreate');
if isempty(pool)
    c = parcluster;                         % get the current cluster object
    numDesiredWorkers = c.NumWorkers;       % get the maximum number of workers
    pool = parpool(numDesiredWorkers);      % start the parallel pool with the maximum number of workers   
end

% Define Simulation Options
sim_id = "berthing_60s_act_iacs";
mkdir(strcat("Results/", sim_id));
sampling_time = 60;                     % seconds
verbose = true;
% final_velocity = -1e-5;                 % -1 cm/s
final_velocity = -5e-6;                 % -5 mm/s
misalignment_type = "oscillating";
state_perturbation_flag = true;
engine_failure_flag = true;
include_actuation = true;

parfor (mc = 1 : MC, pool.NumWorkers)
% for mc = 1 : MC

    try

        fprintf('Starting Simulation n° %2d...\n', mc);

        [RHO_LVLH, M_ctrl_DA, M_ctrl, M_drift, DU, TU, RHOd_LVLH, dist, vel, ...
            renderdata, TCC, Xt_MCI, RHO_MCI, u, u_norms, f_norms, kp_store, ...
            qe0, qe, Tc, Ta, omega_e, omega_e_norms, angle_e, betas, gammas, acc, ...
            deltaState, tspan, tspan_ctrl, Y_ctrl, t0, tf, failure_times, misalignments, ...
            Y_drift, Q_N2C_drift, qe0_drift, qe_drift, Tc_drift, Ta_drift, ...
            omega_e_drift, omega_e_drift_norms, angle_e_drift] = parfmain(sampling_time, include_actuation, final_velocity, verbose, misalignment_type, state_perturbation_flag, engine_failure_flag);

        is_safe = check_min_distance(dist, DU, M_ctrl_DA, M_ctrl, 9.5);

        data(mc).color = cmap(mc, :);
        data(mc).RHO_LVLH = RHO_LVLH;
        data(mc).M_ctrl_DA = M_ctrl_DA;
        data(mc).M_ctrl = M_ctrl;
        data(mc).M_drift = M_drift;
        data(mc).DU = DU;
        data(mc).TU = TU;
        data(mc).RHOd_LVLH = RHOd_LVLH;
        data(mc).dist = dist;
        data(mc).vel = vel;
        data(mc).renderdata = renderdata;
        data(mc).TCC = TCC;
        data(mc).Xt_MCI = Xt_MCI;
        data(mc).RHO_MCI = RHO_MCI;
        data(mc).u = u;
        data(mc).u_norms = u_norms;
        data(mc).f_norms = f_norms;
        data(mc).kp_store = kp_store;
        data(mc).qe0 = qe0;
        data(mc).qe = qe;
        data(mc).Tc = Tc;
        data(mc).Ta = Ta;
        data(mc).omega_e = omega_e;
        data(mc).omega_e_norms = omega_e_norms;
        data(mc).angle_e = angle_e;
        data(mc).betas = betas;
        data(mc).gammas = gammas;
        data(mc).acc = acc;
        data(mc).deltaState = deltaState;
        data(mc).tspan = tspan;
        data(mc).tspan_ctrl = tspan_ctrl;
        data(mc).Y_ctrl = Y_ctrl;
        data(mc).t0 = t0;
        data(mc).tf = tf;
        data(mc).failure_times = failure_times;
        data(mc).misalignments = misalignments;
        data(mc).Y_drift = Y_drift;
        data(mc).Q_N2C_drift = Q_N2C_drift;
        data(mc).qe0_drift = qe0_drift;
        data(mc).qe_drift = qe_drift;
        data(mc).Tc_drift = Tc_drift;
        data(mc).Ta_drift = Ta_drift;
        data(mc).omega_e_drift = omega_e_drift;
        data(mc).omega_e_drift_norms = omega_e_drift_norms;
        data(mc).angle_e_drift = angle_e_drift;

    
        if norm(deltaState(1:3)) <= successful_dist_tol && norm(deltaState(4:6)) <= successful_vel_tol
            data(mc).status = true;
        else
            data(mc).status = false;
        end
    
        if ~is_safe
            data(mc).status = data(mc).status - 0.5;
        end
    
        table(mc, :) = [mc, data(mc).status, deltaState, norm(deltaState(1:3)), norm(deltaState(4:6))];

    catch

        fprintf('Simulation n° %2d was not successful.\n', mc);
        data(mc).status = -1;
        table(mc, :) = [mc, data(mc).status, zeros(1, 8)];

    end


end

% save(strcat("Results/", sim_id, "/", sim_id, ".mat"), 'sim_id', 'sampling_time', 'final_velocity', 'misalignment_type', 'engine_failure_flag', 'state_perturbation_flag', 'include_actuation', 'data', 'table', 'MC');
save(strcat("Results/", sim_id, "/", sim_id, ".mat"));
