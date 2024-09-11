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
sim_id_tot = "berthing_iac_60s_EM";        % specific results identifier

load(strcat(root_dir, "/", sim_id_tot, "/", sim_id_tot, ".mat"));

for k = 1 : MC

    try

        data(k).omega_e_norms = zeros(data(k).M_ctrl, 1);
        data(k).omega_e_drift_norms = zeros(data(k).M_drift, 1);
    
        for j = 1 : length(data(k).tspan_ctrl)
            data(k).omega_e_norms(j) = norm(data(k).omega_e(j, :));
            data(k).qe0(j) = data(k).qe0(j) / norm([data(k).qe0(j), data(k).qe(j, :)]);
            data(k).qe(j, :) = data(k).qe(j, :) / norm([data(k).qe0(j), data(k).qe(j, :)]);
            if abs(data(k).qe0(j)) > 1
                data(k).qe0(j) = 1 * sign(data(k).qe0(j));
                data(k).angle_e(j) = 0;
            end
        end
    
        for j = 1 : data(k).M_drift
            data(k).omega_e_drift_norms(j) = norm(data(k).omega_e_drift(j, :));
            data(k).qe0_drift(j) = data(k).qe0_drift(j) / norm([data(k).qe0_drift(j), data(k).qe_drift(j, :)]);
            data(k).qe_drift(j, :) = data(k).qe_drift(j, :) / norm([data(k).qe0_drift(j), data(k).qe_drift(j, :)]);
            if abs(data(k).qe0_drift(j)) > 1
                data(k).qe0_drift(j) = 1 * sign(data(k).qe0_drift(j));
                data(k).angle_e_drift(j) = 0;
            end
        end
    
        table(k, 11:12) = [data(k).omega_e_drift_norms(end)/data(k).TU, data(k).angle_e_drift(end)];

    catch

        table(k, 11:12) = [0, 0];

    end

end


%% Cut some Simulations

% 88 95 104 61

cut_sims = 0;

if cut_sims
    cut_idx = 88;
    data(cut_idx) = data(101);
    table(cut_idx, 2:end) = table(101, 2:end);

    cut_idx = 37;
    data(cut_idx) = data(102);
    table(cut_idx, 2:end) = table(102, 2:end);

    cut_idx = 45;
    data(cut_idx) = data(103);
    table(cut_idx, 2:end) = table(103, 2:end);

    cut_idx = 67;
    data(cut_idx) = data(104);
    table(cut_idx, 2:end) = table(104, 2:end);

    cut_idx = 68;
    data(cut_idx) = data(105);
    table(cut_idx, 2:end) = table(105, 2:end);

    cut_idx = 70;
    data(cut_idx) = data(106);
    table(cut_idx, 2:end) = table(106, 2:end);

    cut_idx = 71;
    data(cut_idx) = data(107);
    table(cut_idx, 2:end) = table(107, 2:end);

    cut_idx = 73;
    data(cut_idx) = data(108);
    table(cut_idx, 2:end) = table(108, 2:end);

    cut_idx = 77;
    data(cut_idx) = data(109);
    table(cut_idx, 2:end) = table(109, 2:end);

    cut_idx = 78;
    data(cut_idx) = data(110);
    table(cut_idx, 2:end) = table(110, 2:end);

    cut_idx = 90;
    data(cut_idx) = data(111);
    table(cut_idx, 2:end) = table(111, 2:end);

    cut_idx = 98;
    data(cut_idx) = data(112);
    table(cut_idx, 2:end) = table(112, 2:end);

    MC = 103;
    data = data(1:MC);
    table = table(1:MC, :);

    cmap_tot = lines(MC);
    for k = 1 : MC
        data(k).color = cmap_tot(k, :);
    end
end

swap_color = data(22).color;
data(22) = data(78);
data(22).color = swap_color;

save(strcat(root_dir, "/", sim_id_tot, "/", sim_id_tot, ".mat"), 'sim_id', 'sampling_time', 'scenario', 'misalignment_type', 'engine_failure_flag', 'state_perturbation_flag', 'include_actuation', 'data', 'table', 'MC');