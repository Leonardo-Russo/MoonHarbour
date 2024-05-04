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
sim_id = "engine_failure";    % specific results identifier


%% Extract the Results

% Get all files for the simulation
sim_dir = strcat(root_dir, "/", sim_id, "/");
files = dir(fullfile(sim_dir, '*.mat'));
files = files(~strcmp({files.name}, 'data.mat'));
n_sims = length(files);

% Preallocate Data
data = struct('RHO_LVLH', [], 'M_ctrl_DA', [], 'DU', [], 'RHOd_LVLH', [], 'color', [], 'dist', [], 'vel', [], 'successful', [], 'deltaState', [], 'failure_times', [], 'misalignments', []);
data = repmat(data, n_sims, 1);

% Preallocate Table
table = zeros(n_sims, 8);

% Define colormap
cmap = lines(n_sims);

successful_dist_tol = 1e-1;     % 10 cm
successful_vel_tol = 1e-2;      % 1 cm/s
success_color = '#39e617';
failure_color = '#d12828';


for k = 1 : n_sims

    % Extract the simulation index from the filename
    filename = files(k).name;
    index_str = regexp(filename, '\d+', 'match');
    index = str2double(index_str{1});  % Converts the string to a double
    
    try
        
        temp = load(fullfile(files(k).folder, files(k).name), 'RHO_LVLH', 'M_ctrl_DA', 'DU', 'RHOd_LVLH', 'dist', 'vel', 'deltaState', 'failure_times', 'misalignments');

        data(index).color = cmap(index, :);
        data(index).RHO_LVLH = temp.RHO_LVLH;
        data(index).M_ctrl_DA = temp.M_ctrl_DA;
        data(index).DU = temp.DU;
        data(index).RHOd_LVLH = temp.RHOd_LVLH;
        data(index).dist = temp.dist;
        data(index).vel = temp.vel;
    
        if norm(temp.deltaState(1:3)) <= successful_dist_tol && norm(temp.deltaState(4:6)) <= successful_vel_tol
            data(index).status = true;
        else
            data(index).status = false;
        end
    
        table(index, :) = [index, data(index).status, temp.deltaState];

    catch

        fprintf('Simulation n° %2d was not accessible.\n', index);
        data(index).status = -1;
        table(index, 1:2) = [index, data(index).status];

    end
    
end

save(strcat(sim_dir, "data.mat"));


%% Visualize the Simulations

close all

% Create the Results table
results_table = array2table(table, 'VariableNames', {'id', 'status', 'dr (m)', 'dtheta (m)', 'dh (m)', 'dv_r (m/s)', 'dv_theta (m/s)', 'dv_h (m/s)'});
excel_filepath = fullfile(sim_dir, strcat(sim_id, ".xlsx"));
writetable(results_table, excel_filepath);
disp(results_table);
fprintf('Results have been saved to: "%s"\n', excel_filepath);

terminal_traj = figure('name', 'Terminal Chaser Trajectory in LVLH Space', 'WindowState', 'maximized');
% title('Terminal Chaser LVLH Trajectory')
for k = 1 : n_sims

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
savefig(terminal_traj, fullfile(sim_dir, strcat(sim_id, ".fig")));
print(terminal_traj, fullfile(sim_dir, strcat(sim_id, ".png")), '-dpng', '-r1000');          % 1000 DPI

