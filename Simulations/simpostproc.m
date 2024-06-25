%% Simulation Post-Processing

close all
clear 
clc

addpath('../')
addpath('../Library/')
addpath('../Data/')
addpath('../Data/Planets/')
addpath('../Data/Materials/')
addpath('../Data/Ephemeris/')

%% Post-Processing

root_dir = "Results";       % root results folder
sim_id_tot = "berthing_60s_act_new_part1";        % specific results identifier

load(strcat(root_dir, "/", sim_id_tot, "/", sim_id_tot, ".mat"));

for k = 1 : MC

    data(k).omega_e_norms = zeros(data(k).M_ctrl, 1);
    data(k).omega_e_drift_norms = zeros(data(k).M_drift, 1);

    for j = 1 : length(data(k).tspan_ctrl)
        data(k).omega_e_norms(j) = norm(data(k).omega_e(j, :));
        if abs(data(k).qe0(j)) > 1
            data(k).qe0(j) = 1 * sign(data(k).qe0(j));
            data(k).angle_e(j) = 0;
        end
    end

    for j = 1 : data(k).M_drift
        data(k).omega_e_drift_norms(j) = norm(data(k).omega_e_drift(j, :));
        if abs(data(k).qe0_drift(j)) > 1
            data(k).qe0_drift(j) = 1 * sign(data(k).qe0_drift(j));
            data(k).angle_e_drift(j) = 0;
        end
    end

    table(k, 11:12) = [data(k).omega_e_drift_norms(end)/data(k).TU, data(k).angle_e_drift(end)];

end

save(strcat(root_dir, "/", sim_id_tot, "/", sim_id_tot, ".mat"), 'sim_id', 'sampling_time', 'final_velocity', 'misalignment_type', 'engine_failure_flag', 'state_perturbation_flag', 'include_actuation', 'data', 'table', 'MC');