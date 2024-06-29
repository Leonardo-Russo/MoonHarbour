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
sim_id_tot = "berthing_60s_act";        % specific results identifier

% Docking 60s
% 1) 54 is fucked
% 2) 52 is fucked

% Docking 60s act
% 1) -
% 2) -

% Berthing 60s
% 1) 12 is the emergency one

% Berthing 60s act
% 1) -
% 2) 30 is the emergency one

%% Extract the Results

load(strcat(root_dir, "/", sim_id_tot, "/", sim_id_tot, ".mat"));

res = '-r600';
fontsize_labels = 24;
fontsize_axes = 24;
line_width = 2;
fontsize_labels = 12;
fontsize_axes = 12;
line_width = 1.2;

close all

% Create the Results table
results_table = array2table(table, 'VariableNames', {'id', 'status', 'dr (m)', 'dtheta (m)', 'dh (m)', 'dv_r (m/s)', 'dv_theta (m/s)', 'dv_h (m/s)', 'delta_rho (m)', 'delta_rhodot (m/s)', 'omega_ef (rad/s)', 'angle_e (deg)'});
excel_filepath = fullfile(root_dir, sim_id_tot, strcat(sim_id_tot, ".xlsx"));
writetable(results_table, excel_filepath);
disp(results_table);
fprintf('Results have been saved to: "%s"\n', excel_filepath);

fprintf('\nFinal Position Error:\nmean = %.4f mm    std = %.4f mm\n\nFinal Velocity Error:\nmean = %.4f mm/s    std = %.4f mm/s\n\nFinal omega_e:\nmean = %.6fe-6 rad/s  std = %.6fe-6 rad/s\n\nFinal angle_e:\nmean = %.4f deg    std = %.4f deg\n\n', mean(table(:, 9))*1e3, std(table(:, 9))*1e3, mean(table(:, 10))*1e3, std(table(:, 10))*1e3, mean(table(:, 11))*1e6, std(table(:, 11))*1e6, rad2deg(mean(table(:, 12))), rad2deg(std(table(:, 12))))

is_safe = zeros(MC, 1);
for k = 1 : MC

    if data(k).status == -1
        continue;                   % skip failed simulations
    end
    
    min_dist = 9.8;
    is_safe(k) = checkMultipleCrossings(data(k).dist*data(k).DU*1e3, min_dist);
    if is_safe(k) > 0
        fprintf('Warning: Distance < %1.1fm in %.0f\n', min_dist, k);
    end

end
% return
terminal_traj = figure('name', 'Terminal Chaser Trajectory in LVLH Space', 'WindowState', 'maximized');
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

    DrawTrajLVLH3D(data(k).RHO_LVLH(data(k).M_ctrl_DA+1:end, 1:3) * data(k).DU, color, '-', false);

end
[x,y,z]=sphere;
rT = 5e-3;      % km - S/C approximated as a sphere of 5 meter radius
% I = imread('titanium.jpg');
I = imread('black.jpg');
surface(rT*x, rT*y, rT*z, flipud(I), 'FaceColor', 'texturemap', 'EdgeColor', 'none', 'CDataMapping', 'direct')
hold on
radius = 10e-3;
surface(radius*x, radius*y, radius*z, 'FaceColor', '#ffcf82', 'EdgeColor', 'none', 'CDataMapping', 'direct', 'FaceAlpha', 0.2)
view(-55, 15)
set(gca, 'FontSize', fontsize_axes);
savefig(terminal_traj, fullfile(root_dir, sim_id_tot, strcat(sim_id_tot, ".fig")));
print(terminal_traj, fullfile(root_dir, sim_id_tot, strcat(sim_id_tot, ".png")), '-dpng', res);          % 300 DPI


r_comp = figure('Name', 'State Components - r', 'WindowState', 'maximized');
for k = 1 : MC

    if data(k).status == -1
        continue;                   % skip failed simulations
    end

    color = data(k).color;

    plot((data(k).tspan(data(k).M_ctrl_DA+1:end) - data(k).t0)*data(k).TU/3600, data(k).RHO_LVLH(data(k).M_ctrl_DA+1:end, 1)*data(k).DU, 'Color', data(k).color, 'LineStyle','-', 'LineWidth', line_width);
    hold on
    grid on

end
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
ylabel('$\rho_r \ [km]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
xlim([(data(1).tspan(data(1).M_ctrl_DA+1) - data(1).t0)*data(1).TU/3600, (data(1).tspan(end) - data(1).t0)*data(1).TU/3600])
set(gca, 'FontSize', fontsize_axes);
savefig(r_comp, fullfile(root_dir, sim_id_tot, "state_components - r.fig"));
print(r_comp, fullfile(root_dir, sim_id_tot, "state_components - r.png"), '-dpng', res);          % 300 DPI

theta_comp = figure('Name', 'State Components - theta', 'WindowState', 'maximized');
for k = 1 : MC

    if data(k).status == -1
        continue;                   % skip failed simulations
    end

    color = data(k).color;

    plot((data(k).tspan(data(k).M_ctrl_DA+1:end) - data(k).t0)*data(k).TU/3600, data(k).RHO_LVLH(data(k).M_ctrl_DA+1:end, 2)*data(k).DU, 'Color', data(k).color, 'LineStyle','-', 'LineWidth', line_width);
    hold on
    grid on

end
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
ylabel('$\rho_\theta \ [km]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
xlim([(data(1).tspan(data(1).M_ctrl_DA+1) - data(1).t0)*data(1).TU/3600, (data(1).tspan(end) - data(1).t0)*data(1).TU/3600])
set(gca, 'FontSize', fontsize_axes);
savefig(theta_comp, fullfile(root_dir, sim_id_tot, "state_components - theta.fig"));
print(theta_comp, fullfile(root_dir, sim_id_tot, "state_components - theta.png"), '-dpng', res);          % 300 DPI

h_comp = figure('Name', 'State Components - h', 'WindowState', 'maximized');
for k = 1 : MC

    if data(k).status == -1
        continue;                   % skip failed simulations
    end

    color = data(k).color;

    plot((data(k).tspan(data(k).M_ctrl_DA+1:end) - data(k).t0)*data(k).TU/3600, data(k).RHO_LVLH(data(k).M_ctrl_DA+1:end, 3)*data(k).DU, 'Color', data(k).color, 'LineStyle','-', 'LineWidth', line_width);
    hold on
    grid on

end
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
ylabel('$\rho_h \ [km]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
xlim([(data(1).tspan(data(1).M_ctrl_DA+1) - data(1).t0)*data(1).TU/3600, (data(1).tspan(end) - data(1).t0)*data(1).TU/3600])
set(gca, 'FontSize', fontsize_axes);
savefig(h_comp, fullfile(root_dir, sim_id_tot, "state_components - h.fig"));
print(h_comp, fullfile(root_dir, sim_id_tot, "state_components - h.png"), '-dpng', res);          % 300 DPI

rdot_comp = figure('Name', 'State Components - rdot', 'WindowState', 'maximized');
for k = 1 : MC

    if data(k).status == -1
        continue;                   % skip failed simulations
    end

    color = data(k).color;

    plot((data(k).tspan(data(k).M_ctrl_DA+1:end) - data(k).t0)*data(k).TU/3600, data(k).RHO_LVLH(data(k).M_ctrl_DA+1:end, 4)*data(k).DU, 'Color', data(k).color, 'LineStyle','-', 'LineWidth', line_width);
    hold on
    grid on

end
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
ylabel('$\dot{\rho}_r \ [km]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
xlim([(data(1).tspan(data(1).M_ctrl_DA+1) - data(1).t0)*data(1).TU/3600, (data(1).tspan(end) - data(1).t0)*data(1).TU/3600])
set(gca, 'FontSize', fontsize_axes);
savefig(rdot_comp, fullfile(root_dir, sim_id_tot, "state_components - rdot.fig"));
print(rdot_comp, fullfile(root_dir, sim_id_tot, "state_components - rdot.png"), '-dpng', res);          % 300 DPI

thetadot_comp = figure('Name', 'State Components - thetadot', 'WindowState', 'maximized');
for k = 1 : MC

    if data(k).status == -1
        continue;                   % skip failed simulations
    end

    color = data(k).color;

    plot((data(k).tspan(data(k).M_ctrl_DA+1:end) - data(k).t0)*data(k).TU/3600, data(k).RHO_LVLH(data(k).M_ctrl_DA+1:end, 5)*data(k).DU, 'Color', data(k).color, 'LineStyle','-', 'LineWidth', line_width);
    hold on
    grid on

end
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
ylabel('$\dot{\rho}_\theta \ [km]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
xlim([(data(1).tspan(data(1).M_ctrl_DA+1) - data(1).t0)*data(1).TU/3600, (data(1).tspan(end) - data(1).t0)*data(1).TU/3600])
set(gca, 'FontSize', fontsize_axes);
savefig(thetadot_comp, fullfile(root_dir, sim_id_tot, "state_components - thetadot.fig"));
print(thetadot_comp, fullfile(root_dir, sim_id_tot, "state_components - thetadot.png"), '-dpng', res);          % 300 DPI

hdot_comp = figure('Name', 'State Components - hdot', 'WindowState', 'maximized');
for k = 1 : MC

    if data(k).status == -1
        continue;                   % skip failed simulations
    end

    color = data(k).color;

    plot((data(k).tspan(data(k).M_ctrl_DA+1:end) - data(k).t0)*data(k).TU/3600, data(k).RHO_LVLH(data(k).M_ctrl_DA+1:end, 6)*data(k).DU, 'Color', data(k).color, 'LineStyle','-', 'LineWidth', line_width);
    hold on
    grid on

end
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
ylabel('$\dot{\rho}_h \ [km]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
xlim([(data(1).tspan(data(1).M_ctrl_DA+1) - data(1).t0)*data(1).TU/3600, (data(1).tspan(end) - data(1).t0)*data(1).TU/3600])
set(gca, 'FontSize', fontsize_axes);
savefig(hdot_comp, fullfile(root_dir, sim_id_tot, "state_components - hdot.fig"));
print(hdot_comp, fullfile(root_dir, sim_id_tot, "state_components - hdot.png"), '-dpng', res);          % 300 DPI


qe0_fig = figure('Name', 'Realignment - Error Quaternions - qe0', 'WindowState', 'maximized');
for k = 1 : MC

    if data(k).status == -1
        continue;                   % skip failed simulations
    end

    color = data(k).color;

    plot((data(k).tspan(data(k).M_ctrl_DA+data(k).M_ctrl+1:end) - data(k).t0)*data(k).TU/3600, data(k).qe0_drift, 'Color', data(k).color, 'LineStyle','-', 'LineWidth', line_width);
    hold on
    grid on

end
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
ylabel('$q_{e0}$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
xlim([(data(1).tspan(data(1).M_ctrl_DA+data(1).M_ctrl+1) - data(1).t0)*data(1).TU/3600, (data(1).tspan(end) - data(1).t0)*data(1).TU/3600])
set(gca, 'FontSize', fontsize_axes);
savefig(qe0_fig, fullfile(root_dir, sim_id_tot, "realignment - qe0.fig"));
print(qe0_fig, fullfile(root_dir, sim_id_tot, "realignment - qe0.png"), '-dpng', res);


qe1_fig = figure('Name', 'Realignment - Error Quaternions - qe1', 'WindowState', 'maximized');
for k = 1 : MC

    if data(k).status == -1
        continue;                   % skip failed simulations
    end

    color = data(k).color;

    plot((data(k).tspan(data(k).M_ctrl_DA+data(k).M_ctrl+1:end) - data(k).t0)*data(k).TU/3600, data(k).qe_drift(:, 1), 'Color', data(k).color, 'LineStyle','-', 'LineWidth', line_width);
    hold on
    grid on

end
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
ylabel('$q_{e1}$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
xlim([(data(1).tspan(data(1).M_ctrl_DA+data(1).M_ctrl+1) - data(1).t0)*data(1).TU/3600, (data(1).tspan(end) - data(1).t0)*data(1).TU/3600])
set(gca, 'FontSize', fontsize_axes);
savefig(qe1_fig, fullfile(root_dir, sim_id_tot, "realignment - qe1.fig"));
print(qe1_fig, fullfile(root_dir, sim_id_tot, "realignment - qe1.png"), '-dpng', res);


qe2_fig = figure('Name', 'Realignment - Error Quaternions - qe2', 'WindowState', 'maximized');
for k = 1 : MC

    if data(k).status == -1
        continue;                   % skip failed simulations
    end

    color = data(k).color;

    plot((data(k).tspan(data(k).M_ctrl_DA+data(k).M_ctrl+1:end) - data(k).t0)*data(k).TU/3600, data(k).qe_drift(:, 2), 'Color', data(k).color, 'LineStyle','-', 'LineWidth', line_width);
    hold on
    grid on

end
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
ylabel('$q_{e2}$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
xlim([(data(1).tspan(data(1).M_ctrl_DA+data(1).M_ctrl+1) - data(1).t0)*data(1).TU/3600, (data(1).tspan(end) - data(1).t0)*data(1).TU/3600])
set(gca, 'FontSize', fontsize_axes);
savefig(qe2_fig, fullfile(root_dir, sim_id_tot, "realignment - qe2.fig"));
print(qe2_fig, fullfile(root_dir, sim_id_tot, "realignment - qe2.png"), '-dpng', res);


qe3_fig = figure('Name', 'Realignment - Error Quaternions - qe3', 'WindowState', 'maximized');
for k = 1 : MC

    if data(k).status == -1
        continue;                   % skip failed simulations
    end

    color = data(k).color;

    plot((data(k).tspan(data(k).M_ctrl_DA+data(k).M_ctrl+1:end) - data(k).t0)*data(k).TU/3600, data(k).qe_drift(:, 3), 'Color', data(k).color, 'LineStyle','-', 'LineWidth', line_width);
    hold on
    grid on

end
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
ylabel('$q_{e3}$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
xlim([(data(1).tspan(data(1).M_ctrl_DA+data(1).M_ctrl+1) - data(1).t0)*data(1).TU/3600, (data(1).tspan(end) - data(1).t0)*data(1).TU/3600])
set(gca, 'FontSize', fontsize_axes);
savefig(qe3_fig, fullfile(root_dir, sim_id_tot, "realignment - qe3.fig"));
print(qe3_fig, fullfile(root_dir, sim_id_tot, "realignment - qe3.png"), '-dpng', res);


omega_e1_fig = figure('Name', 'Realignment - Error Quaternions - omega_e1', 'WindowState', 'maximized');
for k = 1 : MC

    if data(k).status == -1
        continue;                   % skip failed simulations
    end

    color = data(k).color;

    plot((data(k).tspan(data(k).M_ctrl_DA+data(k).M_ctrl+1:end) - data(k).t0)*data(k).TU/3600, data(k).omega_e_drift(:, 1)/data(k).TU, 'Color', data(k).color, 'LineStyle','-', 'LineWidth', line_width);
    hold on
    grid on

end
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
ylabel('$\omega_{e1} \ [rad/s]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
xlim([(data(1).tspan(data(1).M_ctrl_DA+data(1).M_ctrl+1) - data(1).t0)*data(1).TU/3600, (data(1).tspan(end) - data(1).t0)*data(1).TU/3600])
set(gca, 'FontSize', fontsize_axes);
savefig(omega_e1_fig, fullfile(root_dir, sim_id_tot, "realignment - omega_e1.fig"));
print(omega_e1_fig, fullfile(root_dir, sim_id_tot, "realignment - omega_e1.png"), '-dpng', res);


omega_e2_fig = figure('Name', 'Realignment - Error Quaternions - omega_e2', 'WindowState', 'maximized');
for k = 1 : MC

    if data(k).status == -1
        continue;                   % skip failed simulations
    end

    color = data(k).color;

    plot((data(k).tspan(data(k).M_ctrl_DA+data(k).M_ctrl+1:end) - data(k).t0)*data(k).TU/3600, data(k).omega_e_drift(:, 2)/data(k).TU, 'Color', data(k).color, 'LineStyle','-', 'LineWidth', line_width);
    hold on
    grid on

end
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
ylabel('$\omega_{e2} \ [rad/s]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
xlim([(data(1).tspan(data(1).M_ctrl_DA+data(1).M_ctrl+1) - data(1).t0)*data(1).TU/3600, (data(1).tspan(end) - data(1).t0)*data(1).TU/3600])
set(gca, 'FontSize', fontsize_axes);
savefig(omega_e2_fig, fullfile(root_dir, sim_id_tot, "realignment - omega_e2.fig"));
print(omega_e2_fig, fullfile(root_dir, sim_id_tot, "realignment - omega_e2.png"), '-dpng', res);


omega_e3_fig = figure('Name', 'Realignment - Error Quaternions - omega_e3', 'WindowState', 'maximized');
for k = 1 : MC

    if data(k).status == -1
        continue;                   % skip failed simulations
    end

    color = data(k).color;

    plot((data(k).tspan(data(k).M_ctrl_DA+data(k).M_ctrl+1:end) - data(k).t0)*data(k).TU/3600, data(k).omega_e_drift(:, 3)/data(k).TU, 'Color', data(k).color, 'LineStyle','-', 'LineWidth', line_width);
    hold on
    grid on

end
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
ylabel('$\omega_{e3} \ [rad/s]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
xlim([(data(1).tspan(data(1).M_ctrl_DA+data(1).M_ctrl+1) - data(1).t0)*data(1).TU/3600, (data(1).tspan(end) - data(1).t0)*data(1).TU/3600])
set(gca, 'FontSize', fontsize_axes);
savefig(omega_e3_fig, fullfile(root_dir, sim_id_tot, "realignment - omega_e3.fig"));
print(omega_e3_fig, fullfile(root_dir, sim_id_tot, "realignment - omega_e3.png"), '-dpng', res);


omega_e_norms_fig = figure('Name', 'Realignment - Error Quaternions - omega_e_norms', 'WindowState', 'maximized');
for k = 1 : MC

    if data(k).status == -1
        continue;                   % skip failed simulations
    end

    color = data(k).color;

    plot((data(k).tspan(data(k).M_ctrl_DA+data(k).M_ctrl+1:end) - data(k).t0)*data(k).TU/3600, data(k).omega_e_drift_norms/data(k).TU, 'Color', data(k).color, 'LineStyle','-', 'LineWidth', line_width);
    hold on
    grid on

end
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
ylabel('$|\omega_{e}| \ [rad/s]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
xlim([(data(1).tspan(data(1).M_ctrl_DA+data(1).M_ctrl+1) - data(1).t0)*data(1).TU/3600, (data(1).tspan(end) - data(1).t0)*data(1).TU/3600])
set(gca, 'FontSize', fontsize_axes);
savefig(omega_e_norms_fig, fullfile(root_dir, sim_id_tot, "realignment - omega_e_norms.fig"));
print(omega_e_norms_fig, fullfile(root_dir, sim_id_tot, "realignment - omega_e_norms.png"), '-dpng', res);


angle_e_fig = figure('Name', 'Realignment - Error Quaternions - angle_e', 'WindowState', 'maximized');
for k = 1 : MC

    if data(k).status == -1
        continue;                   % skip failed simulations
    end

    color = data(k).color;

    plot((data(k).tspan(data(k).M_ctrl_DA+data(k).M_ctrl+1:end) - data(k).t0)*data(k).TU/3600, rad2deg(data(k).angle_e_drift), 'Color', data(k).color, 'LineStyle','-', 'LineWidth', line_width);
    hold on
    grid on

end
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
ylabel('$\Phi_{e} \ [deg]$', 'interpreter', 'latex', 'fontsize', fontsize_labels)
xlim([(data(1).tspan(data(1).M_ctrl_DA+data(1).M_ctrl+1) - data(1).t0)*data(1).TU/3600, (data(1).tspan(end) - data(1).t0)*data(1).TU/3600])
set(gca, 'FontSize', fontsize_axes);
savefig(angle_e_fig, fullfile(root_dir, sim_id_tot, "realignment - angle_e.fig"));
print(angle_e_fig, fullfile(root_dir, sim_id_tot, "realignment - angle_e.png"), '-dpng', res);



