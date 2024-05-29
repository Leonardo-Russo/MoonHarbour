%% Builder for MC Simulations - Leonardo Russo

close all
clear
clc

addpath('../')
addpath('../Library/')
addpath('../Data/')
addpath('../Data/Planets/')
addpath('../Data/Materials/')
addpath('../Data/Ephemeris/')

root_dir = "Results";       % root results folder
sim_id = "combined_60s";        % specific results identifier


%% Extract the Results

load(strcat(root_dir, "/", sim_id, "/", sim_id, ".mat"));

close all

% Create the Results table
results_table = array2table(table, 'VariableNames', {'id', 'status', 'dr (m)', 'dtheta (m)', 'dh (m)', 'dv_r (m/s)', 'dv_theta (m/s)', 'dv_h (m/s)', 'delta_rho (m)', 'delta_rhodot (m/s)'});
excel_filepath = fullfile(root_dir, sim_id, strcat(sim_id, ".xlsx"));
writetable(results_table, excel_filepath);
disp(results_table);
fprintf('Results have been saved to: "%s"\n', excel_filepath);

fprintf('\nFinal Position Error:\nmean = %.6f mm    std = %.6f mm\n\nFinal Velocity Error:\nmean = %.6f mm/s    std = %.6f mm/s\n', mean(table(:, 9))*1e3, std(table(:, 9))*1e3, mean(table(:, 10))*1e3, std(table(:, 10))*1e3)

is_safe = zeros(MC, 1);
for k = 1 : MC
    
    min_dist = 9.8;
    is_safe(k) = checkMultipleCrossings(data(k).dist*data(k).DU*1e3, min_dist);
    if is_safe(k) > 0
        fprintf('Warning: Distance < %1.1fm in %.0f\n', min_dist, k);
    end

end

return
terminal_traj = figure('name', 'Terminal Chaser Trajectory in LVLH Space', 'WindowState', 'maximized');
% title('Terminal Chaser LVLH Trajectory')
for k = 1 : MC

    if data(k).status == -1
        continue;                   % skip failed simulations
    end

    % if data(k).status
    %     color = success_color;
    % else
    %     color = failure_color;
    % end

    color = data(k).color;

    DrawTrajLVLH3D(data(k).RHO_LVLH(data(k).M_ctrl_DA:end, 1:3) * data(k).DU, color);

end
view(-65, 15)

% Save the Figure
% savefig(terminal_traj, fullfile(root_dir, sim_dir, strcat(sim_id, ".fig")));
print(terminal_traj, fullfile(root_dir, sim_id, strcat(sim_id, ".png")), '-dpng', '-r300');          % 300 DPI

