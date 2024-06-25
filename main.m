%% Moon Harbour Project - Leonardo Russo

close all
clear
clc
                                                                                
addpath('Library/')
addpath('Data/')
addpath('Data/Planets/')
addpath('Data/Materials/')
addpath('Data/Ephemeris/')


%% Perform the Simulation

sampling_time = 60;                     % seconds
include_actuation = true;
final_velocity = -1e-5;                 % -1 cm/s
% final_velocity = -5e-6;                 % -5 mm/s
verbose = true;
misalignment_type = "oscillating";
state_perturbation_flag = true;
engine_failure_flag = true;
workspace_path = "Data/Utils/main.mat";

parfmain(sampling_time, include_actuation, final_velocity, verbose, misalignment_type, state_perturbation_flag, engine_failure_flag, workspace_path);


%% Visualize the Results
clc
close all

% load('Data/Utils/main.mat');

% load('Results/berthing_60s_act_combined/main.mat');

sim_from_mc = true;

sim_id = "berthing_60s_act_new_part1";
load(strcat("Simulations/Results/", sim_id, "/", sim_id, ".mat"), 'data');
sim_n = 9;

[RHO_LVLH, M_ctrl_DA, M_ctrl, M_drift, DU, TU, RHOd_LVLH, dist, vel, ...
    renderdata, TCC, Xt_MCI, RHO_MCI, u, u_norms, f_norms, kp_store, ...
    qe0, qe, Tc, Ta, omega_e, omega_e_norms, angle_e, betas, gammas, acc, ...
    deltaState, tspan, tspan_ctrl, Y_ctrl, t0, tf, failure_times, misalignments, ...
    Y_drift, Q_N2C_drift, qe0_drift, qe_drift, Tc_drift, Ta_drift, ...
    omega_e_drift, omega_e_drift_norms, angle_e_drift] = workspace_from_data(data(sim_n));

sec2hrs = 1/3600;
runtime = 0;            % I don't know the runtime :(
g0 = 9.80665;           % m/s^2
u_limit = 5e-5*g0;      % m/s^2
MU = 1/DU^2;
opt.include_actuation = true;


opt.saveplots = true;
opt.additional_plots = false;
include_realignment_manoeuvre = 1;
direct_approach_results = 1;
rdv_with_drift = 1;

res = '-r600';

if ~sim_from_mc
    % Show Numerical Results
    fprintf('\nTotal Runtime: %.1f s.\n\n', runtime)
    VarNames = {'r (m)', 'theta (m)', 'h (m)', 'v_r (m/s)', 'v_theta (m/s)', 'v_h (m/s)'};
    DesRhoState = array2table(desState, 'VariableNames', VarNames, 'RowNames', {'Desired State'});
    FinalRhoState = array2table(finalState, 'VariableNames', VarNames, 'RowNames', {'Final State'});
    DeltaRhoState = array2table(deltaState, 'VariableNames', VarNames, 'RowNames', {'Delta State'});
    ResultsTable = [DesRhoState; FinalRhoState; DeltaRhoState];      % concatenate the two tables vertically
    disp(ResultsTable); 
end

% Terminal Chaser State in LVLH
fig = figure('name', 'Terminal Chaser Trajectory in LVLH Space');
C_LVLH_T = DrawTrajLVLH3D(RHO_LVLH(M_ctrl_DA+1:end, 1:3)*DU);
Cd_LVLH_T = DrawTrajLVLH3D(RHOd_LVLH(M_ctrl_DA+1:end, 1:3)*DU, '#6efad2', '-.');
% vp_T = plot3(viapoints_T(:, 1)*DU, viapoints_T(:, 2)*DU, viapoints_T(:, 3)*DU, 'color', 'r', 'linestyle', 'none', 'marker', '.', 'markersize', 15);
% title('Terminal Chaser LVLH Trajectory')
% legend([C_LVLH_T, Cd_LVLH_T, vp_T], {'Chaser Trajectory', 'Reference Trajectory', 'Terminal Via Points'}, 'location', 'best')
% legend([C_LVLH_T, Cd_LVLH_T], {'Chaser Trajectory', 'Reference Trajectory'}, 'location', 'best')
% view(-55, 15)
view(0, 90)
delete(Cd_LVLH_T)
if opt.saveplots
    print(fig, 'Output/Plots/Trajectory Terminal LVLH.png', '-dpng', res);
end


% Chaser LVLH State Components
fig1 = figure('name', 'Chaser LVLH State Component: rho_r');
plot((tspan(1:size(RHOd_LVLH, 1))-t0)*TU*sec2hrs, RHOd_LVLH(:, 1)*DU, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
hold on
plot((tspan-t0)*TU*sec2hrs, RHO_LVLH(:, 1)*DU, 'color', '#4195e8', 'LineWidth', 1.5)
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\rho_r \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('Desired', 'Actual', 'location', 'best')
grid on
if opt.saveplots
    print(fig1, 'Output/Plots/State_LVLH_Component_rho_r.png', '-dpng', res);
end
fig2 = figure('name', 'Chaser LVLH State Component: rho_theta');
plot((tspan(1:size(RHOd_LVLH, 1))-t0)*TU*sec2hrs, RHOd_LVLH(:, 2)*DU, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
hold on
plot((tspan-t0)*TU*sec2hrs, RHO_LVLH(:, 2)*DU, 'color', '#4195e8', 'LineWidth', 1.5)
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\rho_\theta \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('Desired', 'Actual', 'location', 'best')
grid on
if opt.saveplots
    print(fig2, 'Output/Plots/State_LVLH_Component_rho_theta.png', '-dpng', res);
end
fig3 = figure('name', 'Chaser LVLH State Component: rho_h');
plot((tspan(1:size(RHOd_LVLH, 1))-t0)*TU*sec2hrs, RHOd_LVLH(:, 3)*DU, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
hold on
plot((tspan-t0)*TU*sec2hrs, RHO_LVLH(:, 3)*DU, 'color', '#4195e8', 'LineWidth', 1.5)
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\rho_h \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('Desired', 'Actual', 'location', 'best')
grid on
if opt.saveplots
    print(fig3, 'Output/Plots/State_LVLH_Component_rho_h.png', '-dpng', res);
end
fig4 = figure('name', 'Chaser LVLH State Component: dot_rho_r');
plot((tspan(1:size(RHOd_LVLH, 1))-t0)*TU*sec2hrs, RHOd_LVLH(:, 4)*DU/TU, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
hold on
plot((tspan-t0)*TU*sec2hrs, RHO_LVLH(:, 4)*DU/TU, 'color', '#4195e8', 'LineWidth', 1.5)
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\dot{\rho}_r \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('Desired', 'Actual', 'location', 'best')
grid on
if opt.saveplots
    print(fig4, 'Output/Plots/State_LVLH_Component_dot_rho_r.png', '-dpng', res);
end
fig5 = figure('name', 'Chaser LVLH State Component: dot_rho_theta');
plot((tspan(1:size(RHOd_LVLH, 1))-t0)*TU*sec2hrs, RHOd_LVLH(:, 5)*DU/TU, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
hold on
plot((tspan-t0)*TU*sec2hrs, RHO_LVLH(:, 5)*DU/TU, 'color', '#4195e8', 'LineWidth', 1.5)
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\dot{\rho}_\theta \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('Desired', 'Actual', 'location', 'best')
grid on
if opt.saveplots
    print(fig5, 'Output/Plots/State_LVLH_Component_dot_rho_theta.png', '-dpng', res);
end
fig6 = figure('name', 'Chaser LVLH State Component: dot_rho_h');
plot((tspan(1:size(RHOd_LVLH, 1))-t0)*TU*sec2hrs, RHOd_LVLH(:, 6)*DU/TU, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
hold on
plot((tspan-t0)*TU*sec2hrs, RHO_LVLH(:, 6)*DU/TU, 'color', '#4195e8', 'LineWidth', 1.5)
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\dot{\rho}_h \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('Desired', 'Actual', 'location', 'best')
grid on
if opt.saveplots
    print(fig6, 'Output/Plots/State_LVLH_Component_dot_rho_h.png', '-dpng', res);
end


% Distance and Velocity Norms
fig = figure('Name', 'Distance and Velocity Norms');
subplot(1, 2, 1)
plot((tspan-t0)*TU*sec2hrs, dist*DU, 'color', '#4195e8', 'LineWidth', 1.5)
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$|\rho| \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
% title('Relative Distance')
grid on
subplot(1, 2, 2)
plot((tspan-t0)*TU*sec2hrs, vel*DU/TU, 'color', '#4195e8', 'LineWidth', 1.5)
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$|\dot{\rho}| \ [km/s]$', 'interpreter', 'latex', 'fontsize', 12)
% title('Relative Velocity')
grid on
if opt.saveplots
    print(fig, 'Output/Plots/Relative Distance and Velocity.png', '-dpng', res);
end


% Control Norm
fig = figure('name', 'Control Thrust');
p1 = plot((tspan-t0)*TU*sec2hrs, u_norms*1000*DU/TU^2, 'Color', '#4195e8', 'LineWidth', 1.5);
hold on
p2 = plot((tspan-t0)*TU*sec2hrs, u_limit*ones(length(tspan), 1), 'r--', 'LineWidth', 1.2);
p3 = plot((tspan-t0)*TU*sec2hrs, f_norms*1000*DU/TU^2, 'Color', '#93faad', 'LineWidth', 1.5);
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$[m/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
% title('Control Norm')
grid on
legend([p1, p2, p3], '$|u|$', '$u_{max}$', '$f$','Location', 'south', 'Fontsize', 12, 'Interpreter', 'latex');
hold off
zoomPos = [0.5, 0.6, 0.3, 0.25]; % set position of the zoomed plot [left, bottom, width, height]
axes('position', zoomPos);
box on  % adds a box around the new axes
plot((tspan-t0)*TU*sec2hrs, u_norms*1000*DU/TU^2, 'Color', '#4195e8', 'LineWidth', 1.5);
hold on
plot((tspan-t0)*TU*sec2hrs, u_limit*ones(length(tspan), 1), 'r--', 'LineWidth', 1.2);
plot((tspan-t0)*TU*sec2hrs, f_norms*1000*DU/TU^2, 'Color', '#93faad', 'LineWidth', 1.5);
grid on
xlim([(tspan(M_ctrl_DA+1)-t0)*TU*sec2hrs, (tf-t0)*TU*sec2hrs]);   % set the x and y limits for the zoomed plot based on the final part of the data
% ylim([-u_limit(1)/2, u_limit(1)/2]);
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 10)
ylabel('$[m/s^2]$', 'interpreter', 'latex', 'fontsize', 10)
if opt.saveplots
    print(fig, 'Output/Plots/Control Norm.png', '-dpng', res);
end


% Control Components
fig = figure('name', 'Control Thrust Components');
u1 = plot((tspan-t0)*TU*sec2hrs, u(:, 1)*1000*DU/TU^2, 'LineWidth', 1.5);
hold on
u2 = plot((tspan-t0)*TU*sec2hrs, u(:, 2)*1000*DU/TU^2, 'LineWidth', 1.5);
u3 = plot((tspan-t0)*TU*sec2hrs, u(:, 3)*1000*DU/TU^2, 'LineWidth', 1.5);
ulim = plot((tspan-t0)*TU*sec2hrs, u_limit*ones(length(tspan), 1), 'r--', 'LineWidth', 1.2);
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$[m/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
% title('Control Components')
grid on
legend([u1, u2, u3, ulim], '$u_r$', '$u_{\theta}$', '$u_h$', '$u_{max}$','Location', 'southeast', 'Fontsize', 12, 'Interpreter','latex');
hold off
zoomPos = [0.5, 0.6, 0.3, 0.25];    % set position of the zoomed plot [left, bottom, width, height]
axes('position', zoomPos);
box on  % adds a box around the new axes
plot((tspan-t0)*TU*sec2hrs, u(:, 1)*1000*DU/TU^2, 'LineWidth', 1.5);
hold on
plot((tspan-t0)*TU*sec2hrs, u(:, 2)*1000*DU/TU^2, 'LineWidth', 1.5);
plot((tspan-t0)*TU*sec2hrs, u(:, 3)*1000*DU/TU^2, 'LineWidth', 1.5);
plot((tspan-t0)*TU*sec2hrs, u_limit*ones(length(tspan), 1), 'r--', 'LineWidth', 1.2);
grid on
xlim([(tspan(M_ctrl_DA+1)-t0)*TU*sec2hrs, (tf-t0)*TU*sec2hrs]); % set the x and y limits to focus on the final part of the plot
% ylim([-u_limit(1)/2, u_limit(1)/2]);
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 10)
ylabel('$[m/s^2]$', 'interpreter', 'latex', 'fontsize', 10)
if opt.saveplots
    print(fig, 'Output/Plots/Control Components.png', '-dpng', res);
end


if direct_approach_results

    % Direct Approach - Terminal Chaser State in LVLH
    fig = figure('name', 'Direct Approach - Terminal Chaser Trajectory in LVLH Space');
    C_LVLH_T = DrawTrajLVLH3D(RHO_LVLH(1:M_ctrl_DA, 1:3)*DU);
    Cd_LVLH_T = DrawTrajLVLH3D(RHOd_LVLH(1:M_ctrl_DA, 1:3)*DU, '#6efad2', '-.');
    % vp_T = plot3(viapoints_T(:, 1)*DU, viapoints_T(:, 2)*DU, viapoints_T(:, 3)*DU, 'color', 'r', 'linestyle', 'none', 'marker', '.', 'markersize', 15);
    % vp_DA = plot3(viapoints_DA(:, 1)*DU, viapoints_DA(:, 2)*DU, viapoints_DA(:, 3)*DU, 'color', 'b', 'linestyle', 'none', 'marker', '.', 'markersize', 15);
    % legend([C_LVLH_T, Cd_LVLH_T, vp_DA], {'Chaser Trajectory', 'Reference Trajectory', 'Via Points'}, 'location', 'best')
    legend([C_LVLH_T, Cd_LVLH_T], {'Chaser Trajectory', 'Reference Trajectory'}, 'location', 'best')
    % view(-55, 15)
    if opt.saveplots
        print(fig, 'Output/Plots/Direct Approach - Trajectory LVLH.png', '-dpng', res);
    end
    
    
    % Direct Approach - Chaser LVLH State Components
    fig1 = figure('name', 'Direct Approach - Chaser LVLH State Component: rho_r');
    plot((tspan(1:M_ctrl_DA) - t0)*TU*sec2hrs, RHOd_LVLH(1:M_ctrl_DA, 1)*DU, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
    hold on
    plot((tspan(1:M_ctrl_DA) - t0)*TU*sec2hrs, RHO_LVLH(1:M_ctrl_DA, 1)*DU, 'color', '#4195e8', 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$\rho_r \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('Desired', 'Actual', 'location', 'best')
    grid on
    if opt.saveplots
        print(fig1, 'Output/Plots/Direct_Approach_State_LVLH_Component_rho_r.png', '-dpng', res);
    end
    fig2 = figure('name', 'Direct Approach - Chaser LVLH State Component: rho_theta');
    plot((tspan(1:M_ctrl_DA) - t0)*TU*sec2hrs, RHOd_LVLH(1:M_ctrl_DA, 2)*DU, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
    hold on
    plot((tspan(1:M_ctrl_DA) - t0)*TU*sec2hrs, RHO_LVLH(1:M_ctrl_DA, 2)*DU, 'color', '#4195e8', 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$\rho_\theta \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('Desired', 'Actual', 'location', 'best')
    grid on
    if opt.saveplots
        print(fig2, 'Output/Plots/Direct_Approach_State_LVLH_Component_rho_theta.png', '-dpng', res);
    end
    fig3 = figure('name', 'Direct Approach - Chaser LVLH State Component: rho_h');
    plot((tspan(1:M_ctrl_DA) - t0)*TU*sec2hrs, RHOd_LVLH(1:M_ctrl_DA, 3)*DU, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
    hold on
    plot((tspan(1:M_ctrl_DA) - t0)*TU*sec2hrs, RHO_LVLH(1:M_ctrl_DA, 3)*DU, 'color', '#4195e8', 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$\rho_h \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('Desired', 'Actual', 'location', 'best')
    grid on
    if opt.saveplots
        print(fig3, 'Output/Plots/Direct_Approach_State_LVLH_Component_rho_h.png', '-dpng', res);
    end
    fig4 = figure('name', 'Direct Approach - Chaser LVLH State Component: dot_rho_r');
    plot((tspan(1:M_ctrl_DA) - t0)*TU*sec2hrs, RHOd_LVLH(1:M_ctrl_DA, 4)*DU/TU, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
    hold on
    plot((tspan(1:M_ctrl_DA) - t0)*TU*sec2hrs, RHO_LVLH(1:M_ctrl_DA, 4)*DU/TU, 'color', '#4195e8', 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$\dot{\rho}_r \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('Desired', 'Actual', 'location', 'best')
    grid on
    if opt.saveplots
        print(fig4, 'Output/Plots/Direct_Approach_State_LVLH_Component_dot_rho_r.png', '-dpng', res);
    end
    fig5 = figure('name', 'Direct Approach - Chaser LVLH State Component: dot_rho_theta');
    plot((tspan(1:M_ctrl_DA) - t0)*TU*sec2hrs, RHOd_LVLH(1:M_ctrl_DA, 5)*DU/TU, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
    hold on
    plot((tspan(1:M_ctrl_DA) - t0)*TU*sec2hrs, RHO_LVLH(1:M_ctrl_DA, 5)*DU/TU, 'color', '#4195e8', 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$\dot{\rho}_\theta \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('Desired', 'Actual', 'location', 'best')
    grid on
    if opt.saveplots
        print(fig5, 'Output/Plots/Direct_Approach_State_LVLH_Component_dot_rho_theta.png', '-dpng', res);
    end
    fig6 = figure('name', 'Direct Approach - Chaser LVLH State Component: dot_rho_h');
    plot((tspan(1:M_ctrl_DA) - t0)*TU*sec2hrs, RHOd_LVLH(1:M_ctrl_DA, 6)*DU/TU, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
    hold on
    plot((tspan(1:M_ctrl_DA) - t0)*TU*sec2hrs, RHO_LVLH(1:M_ctrl_DA, 6)*DU/TU, 'color', '#4195e8', 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$\dot{\rho}_h \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('Desired', 'Actual', 'location', 'best')
    grid on
    if opt.saveplots
        print(fig6, 'Output/Plots/Direct_Approach_State_LVLH_Component_dot_rho_h.png', '-dpng', res);
    end

    
    % Direct Approach - Distance and Velocity Norms
    fig = figure('Name', 'Direct Approach - Distance and Velocity Norms');
    subplot(1, 2, 1)
    plot((tspan(1:M_ctrl_DA) - t0)*TU*sec2hrs, dist(1:M_ctrl_DA)*DU, 'color', '#4195e8', 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$|\rho| \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
    % title('Relative Distance')
    grid on
    subplot(1, 2, 2)
    plot((tspan(1:M_ctrl_DA) - t0)*TU*sec2hrs, vel(1:M_ctrl_DA)*DU/TU, 'color', '#4195e8', 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$|\dot{\rho}| \ [km/s]$', 'interpreter', 'latex', 'fontsize', 12)
    % title('Relative Velocity')
    grid on
    if opt.saveplots
        print(fig, 'Output/Plots/Direct Approach - Relative Distance and Velocity.png', '-dpng', res);
    end
    
    
    % Direct Approach - Control Norm
    fig = figure('name', 'Direct Approach - Control Thrust');
    p1 = plot((tspan(1:M_ctrl_DA) - t0)*TU*sec2hrs, u_norms(1:M_ctrl_DA)*1000*DU/TU^2, 'Color', '#4195e8', 'LineWidth', 1.5);
    hold on
    p2 = plot((tspan(1:M_ctrl_DA) - t0)*TU*sec2hrs, u_limit*ones(length(tspan(1:M_ctrl_DA)), 1), 'r--', 'LineWidth', 1.2);
    p3 = plot((tspan(1:M_ctrl_DA) - t0)*TU*sec2hrs, f_norms(1:M_ctrl_DA)*1000*DU/TU^2, 'Color', '#93faad', 'LineWidth', 1.5);
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$[m/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
    % title('Control Norm')
    grid on
    legend([p1, p2, p3], '$|u|$', '$u_{max}$', '$f$','Location', 'south', 'Fontsize', 12, 'Interpreter', 'latex');
    % xlim([(tspan(M_ctrl_DA+1)-t0)*TU*sec2hrs, (tf-t0)*TU*sec2hrs]);
    if opt.saveplots
        print(fig, 'Output/Plots/Direct Approach - Control Norm.png', '-dpng', res);
    end
    
    
    % Direct Approach - Control Components
    fig = figure('name', 'Direct Approach - Control Thrust Components');
    u1 = plot((tspan(1:M_ctrl_DA) - t0)*TU*sec2hrs, u(1:M_ctrl_DA, 1)*1000*DU/TU^2, 'LineWidth', 1.5);
    hold on
    u2 = plot((tspan(1:M_ctrl_DA) - t0)*TU*sec2hrs, u(1:M_ctrl_DA, 2)*1000*DU/TU^2, 'LineWidth', 1.5);
    u3 = plot((tspan(1:M_ctrl_DA) - t0)*TU*sec2hrs, u(1:M_ctrl_DA, 3)*1000*DU/TU^2, 'LineWidth', 1.5);
    ulim = plot((tspan(1:M_ctrl_DA) - t0)*TU*sec2hrs, u_limit*ones(length(tspan(1:M_ctrl_DA)), 1), 'r--', 'LineWidth', 1.2);
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$[m/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
    % title('Control Components')
    grid on
    legend([u1, u2, u3, ulim], '$u_r$', '$u_{\theta}$', '$u_h$', '$u_{max}$','Location', 'southeast', 'Fontsize', 12, 'Interpreter','latex');
    if opt.saveplots
        print(fig, 'Output/Plots/Direct Approach - Control Components.png', '-dpng', res);
    end

end

if ~sim_from_mc
    % Body and Commanded Attitude Quaternions
    fig = figure('name', "Body and Commanded Attitude");
    subplot(1, 2, 1)
    plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 14), 'LineWidth', 1.5)
    hold on
    plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 15), 'LineWidth', 1.5)
    plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 16), 'LineWidth', 1.5)
    plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 17), 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$q_i$', 'interpreter', 'latex', 'fontsize', 12)
    legend('q_{b0}', 'q_{b1}', 'q_{b2}', 'q_{b3}', 'fontsize', 10, 'location', 'best')
    grid on
    subplot(1, 2, 2)
    plot((tspan_ctrl-t0)*TU*sec2hrs, Q_N2C_AOCS_stack(:, 1), 'LineWidth', 1.5)
    hold on
    plot((tspan_ctrl-t0)*TU*sec2hrs, Q_N2C_AOCS_stack(:, 2), 'LineWidth', 1.5)
    plot((tspan_ctrl-t0)*TU*sec2hrs, Q_N2C_AOCS_stack(:, 3), 'LineWidth', 1.5)
    plot((tspan_ctrl-t0)*TU*sec2hrs, Q_N2C_AOCS_stack(:, 4), 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$q_i$', 'interpreter', 'latex', 'fontsize', 12)
    legend('q_{c0}', 'q_{c1}', 'q_{c2}', 'q_{c3}', 'fontsize', 10, 'location', 'best')
    grid on  
    if opt.saveplots
        print(fig, 'Output/Plots/Body and Commanded Quaternions.png', '-dpng', res);
    end
end

% Error Quaternions
fig = figure('name', "Error Quaternions");
plot((tspan_ctrl-t0)*TU*sec2hrs, qe0, 'LineWidth', 1.5)
hold on
plot((tspan_ctrl-t0)*TU*sec2hrs, qe(:, 1), 'LineWidth', 1.5)
plot((tspan_ctrl-t0)*TU*sec2hrs, qe(:, 2), 'LineWidth', 1.5)
plot((tspan_ctrl-t0)*TU*sec2hrs, qe(:, 3), 'LineWidth', 1.5)
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$q_i$', 'interpreter', 'latex', 'fontsize', 12)
legend('q_{e0}', 'q_{e1}', 'q_{e2}', 'q_{e3}', 'fontsize', 10, 'location', 'best')
grid on
if opt.saveplots
    print(fig, 'Output/Plots/Error Quaternions.png', '-dpng', res);
end


% Error Angular Velocity
fig = figure('name', "Error Angular Velocity");
plot((tspan_ctrl-t0)*TU*sec2hrs, omega_e(:, 1)/TU, 'LineWidth', 1.5)
hold on
plot((tspan_ctrl-t0)*TU*sec2hrs, omega_e(:, 2)/TU, 'LineWidth', 1.5)
plot((tspan_ctrl-t0)*TU*sec2hrs, omega_e(:, 3)/TU, 'LineWidth', 1.5)
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\omega_{ei} \, [rad/s]$', 'interpreter', 'latex', 'fontsize', 12)
legend('\omega_{e1}', '\omega_{e2}', '\omega_{e3}', 'fontsize', 10, 'location', 'best')
grid on
if opt.saveplots
    print(fig, 'Output/Plots/Error Angular Velocity.png', '-dpng', res);
end

% Error Angle
fig = figure('name', "Error Angle");
plot((tspan_ctrl-t0)*TU*sec2hrs, rad2deg(angle_e), 'LineWidth', 1.5)
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\Phi_e \, [deg]$', 'interpreter', 'latex', 'fontsize', 12)
grid on
if opt.saveplots
    print(fig, 'Output/Plots/Error Angle.png', '-dpng', res);
end


% Body and Wheels Angular Velocities
fig = figure('name', "Body and Wheels Angular Velocities");
subplot(1, 2, 1)
plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 18)/TU, 'LineWidth', 1.5)
hold on
plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 19)/TU, 'LineWidth', 1.5)
plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 20)/TU, 'LineWidth', 1.5)
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\omega_i \, [rad/s]$', 'interpreter', 'latex', 'fontsize', 12)
legend('\omega_{1}', '\omega_{2}', '\omega_{3}', 'fontsize', 10, 'location', 'best')
grid on
subplot(1, 2, 2)
plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 21)/TU, 'LineWidth', 1.5)
hold on
plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 22)/TU, 'LineWidth', 1.5)
plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 23)/TU, 'LineWidth', 1.5)
plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 24)/TU, 'LineWidth', 1.5)
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\omega_{si} \, [rad/s]$', 'interpreter', 'latex', 'fontsize', 12)
legend('\omega_{s1}', '\omega_{s2}', '\omega_{s3}', '\omega_{s4}', 'fontsize', 10, 'location', 'best')
grid on
if opt.saveplots
    print(fig, 'Output/Plots/Body and Wheels Angular Velocities.png', '-dpng', res);
end


% Commanded Torque
fig = figure('name', 'Commanded Torque');
plot((tspan_ctrl-t0)*TU*sec2hrs, Tc(:, 1)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
hold on
plot((tspan_ctrl-t0)*TU*sec2hrs, Tc(:, 2)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
plot((tspan_ctrl-t0)*TU*sec2hrs, Tc(:, 3)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$T_{c,i} \ [Nm]$', 'interpreter', 'latex', 'fontsize', 12)
legend('T_{c1}', 'T_{c2}', 'T_{c3}', 'fontsize', 10, 'location', 'best')
grid on
if opt.saveplots
    print(fig, 'Output/Plots/Commanded Torque.png', '-dpng', res);
end


if opt.include_actuation
    % Actual Torque
    fig = figure('name', 'Actual Torque');
    plot((tspan_ctrl-t0)*TU*sec2hrs, Ta(:, 1)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
    hold on
    plot((tspan_ctrl-t0)*TU*sec2hrs, Ta(:, 2)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
    plot((tspan_ctrl-t0)*TU*sec2hrs, Ta(:, 3)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$T_{a,i} \ [Nm]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('T_{a1}', 'T_{a2}', 'T_{a3}', 'fontsize', 10, 'location', 'best')
    grid on
    if opt.saveplots
        print(fig, 'Output/Plots/Actual Torque.png', '-dpng', res);
    end
    
    
    % Actual vs Commanded Torque
    fig = figure('name', 'Actual vs Commanded Torque');
    subplot(1, 3, 1)
    plot((tspan_ctrl-t0)*TU*sec2hrs, Ta(:, 1)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
    hold on
    plot((tspan_ctrl-t0)*TU*sec2hrs, Tc(:, 1)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$T_{1} \ [Nm]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('T_{a1}', 'T_{c1}', 'fontsize', 10, 'location', 'best')
    grid on
    subplot(1, 3, 2)
    plot((tspan_ctrl-t0)*TU*sec2hrs, Ta(:, 2)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
    hold on
    plot((tspan_ctrl-t0)*TU*sec2hrs, Tc(:, 2)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$T_{2} \ [Nm]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('T_{a2}', 'T_{c2}', 'fontsize', 10, 'location', 'best')
    grid on
    subplot(1, 3, 3)
    plot((tspan_ctrl-t0)*TU*sec2hrs, Ta(:, 3)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
    hold on
    plot((tspan_ctrl-t0)*TU*sec2hrs, Tc(:, 3)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$T_{3} \ [Nm]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('T_{a3}', 'T_{c3}', 'fontsize', 10, 'location', 'best')
    grid on
    if opt.saveplots
        print(fig, 'Output/Plots/Actual vs Commanded Torque.png', '-dpng', res);
    end
end


% Gammas and Betas
if misalignment_type ~= "null"
    fig = figure('name', "Misalignment Angles - Beta");
    plot((tspan_ctrl-t0)*TU*sec2hrs, rad2deg(betas), 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$\beta \ [deg]$', 'interpreter', 'latex', 'fontsize', 12)
    grid on
    if opt.saveplots
        print(fig, 'Output/Plots/Beta Angles.png', '-dpng', res);
    end

    fig = figure('name', "Misalignment Angles - Gamma");
    plot((tspan_ctrl-t0)*TU*sec2hrs, rad2deg(gammas), 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$\gamma \ [deg]$', 'interpreter', 'latex', 'fontsize', 12)
    grid on   
    if opt.saveplots
        print(fig, 'Output/Plots/Gamma Angles.png', '-dpng', res);
    end
end


if include_realignment_manoeuvre

    % Realignment - Body and Commanded Attitude Quaternions
    fig = figure('name', "Realignment - Body and Commanded Attitude");
    subplot(1, 2, 1)
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Y_drift(:, 14), 'LineWidth', 1.5)
    hold on
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Y_drift(:, 15), 'LineWidth', 1.5)
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Y_drift(:, 16), 'LineWidth', 1.5)
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Y_drift(:, 17), 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$q_i$', 'interpreter', 'latex', 'fontsize', 12)
    legend('q_{b0}', 'q_{b1}', 'q_{b2}', 'q_{b3}', 'fontsize', 10, 'location', 'best')
    grid on
    subplot(1, 2, 2)
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Q_N2C_drift(:, 1), 'LineWidth', 1.5)
    hold on
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Q_N2C_drift(:, 2), 'LineWidth', 1.5)
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Q_N2C_drift(:, 3), 'LineWidth', 1.5)
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Q_N2C_drift(:, 4), 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$q_i$', 'interpreter', 'latex', 'fontsize', 12)
    legend('q_{c0}', 'q_{c1}', 'q_{c2}', 'q_{c3}', 'fontsize', 10, 'location', 'best')
    grid on  
    if opt.saveplots
        print(fig, 'Output/Plots/Realignment - Body and Commanded Quaternions.png', '-dpng', res);
    end
    
    
    % Realignment - Error Quaternions
    fig = figure('name', "Realignment - Error Quaternions");
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, qe0_drift, 'LineWidth', 1.5)
    hold on
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, qe_drift(:, 1), 'LineWidth', 1.5)
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, qe_drift(:, 2), 'LineWidth', 1.5)
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, qe_drift(:, 3), 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$q_i$', 'interpreter', 'latex', 'fontsize', 12)
    legend('q_{e0}', 'q_{e1}', 'q_{e2}', 'q_{e3}', 'fontsize', 10, 'location', 'best')
    grid on
    if opt.saveplots
        print(fig, 'Output/Plots/Realignment - Error Quaternions.png', '-dpng', res);
    end

    % Realignment - Error Angular Velocity
    fig = figure('name', "Realignment - Error Angular Velocity");
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, omega_e_drift(:, 1)/TU, 'LineWidth', 1.5)
    hold on
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, omega_e_drift(:, 2)/TU, 'LineWidth', 1.5)
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, omega_e_drift(:, 3)/TU, 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$\omega_{ei} \, [rad/s]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('\omega_{e1}', '\omega_{e2}', '\omega_{e3}', 'fontsize', 10, 'location', 'best')
    grid on
    if opt.saveplots
        print(fig, 'Output/Plots/Realignment - Error Angular Velocity.png', '-dpng', res);
    end

    % Realignment - Error Angle
    fig = figure('name', "Realignment - Error Angle");
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, rad2deg(angle_e_drift), 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$\Phi_e \, [deg]$', 'interpreter', 'latex', 'fontsize', 12)
    grid on
    if opt.saveplots
        print(fig, 'Output/Plots/Realignment - Error Angle.png', '-dpng', res);
    end
    
    % Realignment - Body and Wheels Angular Velocities
    fig = figure('name', "Realignment - Body and Wheels Angular Velocities");
    subplot(1, 2, 1)
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Y_drift(:, 18)/TU, 'LineWidth', 1.5)
    hold on
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Y_drift(:, 19)/TU, 'LineWidth', 1.5)
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Y_drift(:, 20)/TU, 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$\omega_i \, [rad/s]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('\omega_{1}', '\omega_{2}', '\omega_{3}', 'fontsize', 10, 'location', 'best')
    grid on
    subplot(1, 2, 2)
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Y_drift(:, 21)/TU, 'LineWidth', 1.5)
    hold on
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Y_drift(:, 22)/TU, 'LineWidth', 1.5)
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Y_drift(:, 23)/TU, 'LineWidth', 1.5)
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Y_drift(:, 24)/TU, 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$\omega_{si} \, [rad/s]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('\omega_{s1}', '\omega_{s2}', '\omega_{s3}', '\omega_{s4}', 'fontsize', 10, 'location', 'best')
    grid on
    if opt.saveplots
        print(fig, 'Output/Plots/Realignment - Body and Wheels Angular Velocities.png', '-dpng', res);
    end
    
    
    % Realignment - Commanded Torque
    fig = figure('name', 'Realignment - Commanded Torque');
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Tc_drift(:, 1)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
    hold on
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Tc_drift(:, 2)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
    plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Tc_drift(:, 3)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$T_{c,i} \ [Nm]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('T_{c1}', 'T_{c2}', 'T_{c3}', 'fontsize', 10, 'location', 'best')
    grid on
    if opt.saveplots
        print(fig, 'Output/Plots/Realignment - Commanded Torque.png', '-dpng', res);
    end

    
    if opt.include_actuation
        % Realignment - Actual Torque
        fig = figure('name', 'Realignment - Actual Torque');
        plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Ta_drift(:, 1)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
        hold on
        plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Ta_drift(:, 2)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
        plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Ta_drift(:, 3)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
        xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$T_{a,i} \ [Nm]$', 'interpreter', 'latex', 'fontsize', 12)
        legend('T_{a1}', 'T_{a2}', 'T_{a3}', 'fontsize', 10, 'location', 'best')
        grid on
        if opt.saveplots
            print(fig, 'Output/Plots/Realignment - Actual Torque.png', '-dpng', res);
        end
        
        
        % Realignment - Actual vs Commanded Torque
        fig = figure('name', 'Realignment - Actual vs Commanded Torque');
        subplot(1, 3, 1)
        plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Ta_drift(:, 1)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
        hold on
        plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Tc_drift(:, 1)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
        xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$T_{1} \ [Nm]$', 'interpreter', 'latex', 'fontsize', 12)
        legend('T_{a1}', 'T_{c1}', 'fontsize', 10, 'location', 'best')
        grid on
        subplot(1, 3, 2)
        plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Ta_drift(:, 2)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
        hold on
        plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Tc_drift(:, 2)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
        xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$T_{2} \ [Nm]$', 'interpreter', 'latex', 'fontsize', 12)
        legend('T_{a2}', 'T_{c2}', 'fontsize', 10, 'location', 'best')
        grid on
        subplot(1, 3, 3)
        plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Ta_drift(:, 3)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
        hold on
        plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Tc_drift(:, 3)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
        xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$T_{3} \ [Nm]$', 'interpreter', 'latex', 'fontsize', 12)
        legend('T_{a3}', 'T_{c3}', 'fontsize', 10, 'location', 'best')
        grid on
        if opt.saveplots
            print(fig, 'Output/Plots/Realignment - Actual vs Commanded Torque.png', '-dpng', res);
        end
    end

    if rdv_with_drift

        % % Total - Body and Commanded Attitude Quaternions
        % fig = figure('name', "Total - Body and Commanded Attitude");
        % subplot(1, 2, 1)
        % plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Y_drift(:, 14), 'LineWidth', 1.5)
        % hold on
        % plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Y_drift(:, 15), 'LineWidth', 1.5)
        % plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Y_drift(:, 16), 'LineWidth', 1.5)
        % plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Y_drift(:, 17), 'LineWidth', 1.5)
        % xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
        % ylabel('$q_i$', 'interpreter', 'latex', 'fontsize', 12)
        % legend('q_{b0}', 'q_{b1}', 'q_{b2}', 'q_{b3}', 'fontsize', 10, 'location', 'best')
        % grid on
        % subplot(1, 2, 2)
        % plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Q_N2C_drift(:, 1), 'LineWidth', 1.5)
        % hold on
        % plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Q_N2C_drift(:, 2), 'LineWidth', 1.5)
        % plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Q_N2C_drift(:, 3), 'LineWidth', 1.5)
        % plot((tspan(M_ctrl_DA+M_ctrl+1:end)-t0)*TU*sec2hrs, Q_N2C_drift(:, 4), 'LineWidth', 1.5)
        % xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
        % ylabel('$q_i$', 'interpreter', 'latex', 'fontsize', 12)
        % legend('q_{c0}', 'q_{c1}', 'q_{c2}', 'q_{c3}', 'fontsize', 10, 'location', 'best')
        % grid on  
        % if opt.saveplots
        %     print(fig, 'Output/Plots/Total - Body and Commanded Quaternions.png', '-dpng', res);
        % end
        
        
        % Total - Error Quaternions
        fig = figure('name', "Total - Error Quaternions");
        plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [qe0; qe0_drift], 'LineWidth', 1.5)
        hold on
        plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [qe(:, 1); qe_drift(:, 1)], 'LineWidth', 1.5)
        plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [qe(:, 2); qe_drift(:, 2)], 'LineWidth', 1.5)
        plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [qe(:, 3); qe_drift(:, 3)], 'LineWidth', 1.5)
        xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$q_i$', 'interpreter', 'latex', 'fontsize', 12)
        legend('q_{e0}', 'q_{e1}', 'q_{e2}', 'q_{e3}', 'fontsize', 10, 'location', 'best')
        grid on
        if opt.saveplots
            print(fig, 'Output/Plots/Total - Error Quaternions.png', '-dpng', res);
        end
    
        % Total - Error Angular Velocity
        fig = figure('name', "Total - Error Angular Velocity");
        plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [omega_e(:, 1); omega_e_drift(:, 1)]/TU, 'LineWidth', 1.5)
        hold on
        plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [omega_e(:, 2); omega_e_drift(:, 2)]/TU, 'LineWidth', 1.5)
        plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [omega_e(:, 3); omega_e_drift(:, 3)]/TU, 'LineWidth', 1.5)
        xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$\omega_{ei} \, [rad/s]$', 'interpreter', 'latex', 'fontsize', 12)
        legend('\omega_{e1}', '\omega_{e2}', '\omega_{e3}', 'fontsize', 10, 'location', 'best')
        grid on
        if opt.saveplots
            print(fig, 'Output/Plots/Total - Error Angular Velocity.png', '-dpng', res);
        end
    
        % Total - Error Angle
        fig = figure('name', "Total - Error Angle");
        plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, rad2deg([angle_e; angle_e_drift]), 'LineWidth', 1.5)
        xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$\Phi_e \, [deg]$', 'interpreter', 'latex', 'fontsize', 12)
        grid on
        if opt.saveplots
            print(fig, 'Output/Plots/Total - Error Angle.png', '-dpng', res);
        end
        
        % Total - Body and Wheels Angular Velocities
        fig = figure('name', "Total - Body and Wheels Angular Velocities");
        subplot(1, 2, 1)
        plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [Y_ctrl(:, 18); Y_drift(:, 18)]/TU, 'LineWidth', 1.5)
        hold on
        plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [Y_ctrl(:, 19); Y_drift(:, 19)]/TU, 'LineWidth', 1.5)
        plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [Y_ctrl(:, 20); Y_drift(:, 20)]/TU, 'LineWidth', 1.5)
        xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$\omega_i \, [rad/s]$', 'interpreter', 'latex', 'fontsize', 12)
        legend('\omega_{1}', '\omega_{2}', '\omega_{3}', 'fontsize', 10, 'location', 'best')
        grid on
        subplot(1, 2, 2)
        plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [Y_ctrl(:, 21); Y_drift(:, 21)]/TU, 'LineWidth', 1.5)
        hold on
        plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [Y_ctrl(:, 22); Y_drift(:, 22)]/TU, 'LineWidth', 1.5)
        plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [Y_ctrl(:, 23); Y_drift(:, 23)]/TU, 'LineWidth', 1.5)
        plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [Y_ctrl(:, 24); Y_drift(:, 24)]/TU, 'LineWidth', 1.5)
        xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$\omega_{si} \, [rad/s]$', 'interpreter', 'latex', 'fontsize', 12)
        legend('\omega_{s1}', '\omega_{s2}', '\omega_{s3}', '\omega_{s4}', 'fontsize', 10, 'location', 'best')
        grid on
        if opt.saveplots
            print(fig, 'Output/Plots/Total - Body and Wheels Angular Velocities.png', '-dpng', res);
        end
        
        
        % Total - Commanded Torque
        fig = figure('name', 'Total - Commanded Torque');
        plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [Tc(:, 1); Tc_drift(:, 1)]*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
        hold on
        plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [Tc(:, 2); Tc_drift(:, 2)]*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
        plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [Tc(:, 3); Tc_drift(:, 3)]*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
        xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$T_{c,i} \ [Nm]$', 'interpreter', 'latex', 'fontsize', 12)
        legend('T_{c1}', 'T_{c2}', 'T_{c3}', 'fontsize', 10, 'location', 'best')
        grid on
        if opt.saveplots
            print(fig, 'Output/Plots/Total - Commanded Torque.png', '-dpng', res);
        end
    
        
        if opt.include_actuation
            % Total - Actual Torque
            fig = figure('name', 'Total - Actual Torque');
            plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [Ta(:, 1); Ta_drift(:, 1)]*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
            hold on
            plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [Ta(:, 2); Ta_drift(:, 2)]*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
            plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [Ta(:, 3); Ta_drift(:, 3)]*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
            xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
            ylabel('$T_{a,i} \ [Nm]$', 'interpreter', 'latex', 'fontsize', 12)
            legend('T_{a1}', 'T_{a2}', 'T_{a3}', 'fontsize', 10, 'location', 'best')
            grid on
            if opt.saveplots
                print(fig, 'Output/Plots/Total - Actual Torque.png', '-dpng', res);
            end
            
            
            % % Total - Actual vs Commanded Torque
            % fig = figure('name', 'Total - Actual vs Commanded Torque');
            % subplot(1, 3, 1)
            % plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [Ta(:, 1); Ta_drift(:, 1)]*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
            % hold on
            % plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [Tc(:, 1); Tc_drift(:, 1)]*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
            % xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
            % ylabel('$T_{1} \ [Nm]$', 'interpreter', 'latex', 'fontsize', 12)
            % legend('T_{a1}', 'T_{c1}', 'fontsize', 10, 'location', 'best')
            % grid on
            % subplot(1, 3, 2)
            % plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [Ta(:, 2); Ta_drift(:, 2)]*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
            % hold on
            % plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [Tc(:, 2); Tc_drift(:, 2)]*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
            % xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
            % ylabel('$T_{2} \ [Nm]$', 'interpreter', 'latex', 'fontsize', 12)
            % legend('T_{a2}', 'T_{c2}', 'fontsize', 10, 'location', 'best')
            % grid on
            % subplot(1, 3, 3)
            % plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [Ta(:, 3); Ta_drift(:, 3)]*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
            % hold on
            % plot((tspan(M_ctrl_DA+1:end)-t0)*TU*sec2hrs, [Tc(:, 3); Tc_drift(:, 3)]*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
            % xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
            % ylabel('$T_{3} \ [Nm]$', 'interpreter', 'latex', 'fontsize', 12)
            % legend('T_{a3}', 'T_{c3}', 'fontsize', 10, 'location', 'best')
            % grid on
            % if opt.saveplots
            %     print(fig, 'Output/Plots/Total - Actual vs Commanded Torque.png', '-dpng', res);
            % end
        end

    end

end


if opt.additional_plots

    % Chaser State in LVLH
    fig = figure('name', 'Chaser Trajectory in LVLH Space');
    C_LVLH = DrawTrajLVLH3D(RHO_LVLH(:, 1:3)*DU);
    Cd_LVLH = DrawTrajLVLH3D(RHOd_LVLH(:, 1:3)*DU, '#6efad2', '-.');
    % vp_T = plot3(viapoints_T(:, 1)*DU, viapoints_T(:, 2)*DU, viapoints_T(:, 3)*DU, 'color', 'r', 'linestyle', 'none', 'marker', '.', 'markersize', 15);
    % vp_DA = plot3(viapoints_DA(:, 1)*DU, viapoints_DA(:, 2)*DU, viapoints_DA(:, 3)*DU, 'color', 'b', 'linestyle', 'none', 'marker', '.', 'markersize', 15);
    % title('Chaser LVLH Trajectory')
    % legend([C_LVLH, Cd_LVLH, vp_T, vp_DA], {'Chaser Trajectory', 'Reference Trajectory', 'Terminal Via Points', 'Direct Approach Via Points'}, 'location', 'best')
    % legend([C_LVLH, Cd_LVLH, vp_DA], {'Chaser Trajectory', 'Reference Trajectory', 'Direct Approach Via Points'}, 'location', 'best')
    legend([C_LVLH, Cd_LVLH], {'Chaser Trajectory', 'Reference Trajectory'}, 'location', 'best')
    if opt.saveplots
        print(fig, 'Output/Plots/Trajectory LVLH.png', '-dpng', res);
    end

    % Target, Chaser and Reference Chaser Trajectories in MCI
    fig = figure('name', 'Trajectory in MCI Space');
    T = DrawTrajMCI3D(Xt_MCI(:, 1:3)*DU, '#d1d1d1', '-.');
    C = DrawTrajMCI3D(Xc_MCI(:, 1:3)*DU);
    legend([T, C], {'Target Trajectory', 'Chaser Trajectory'}, 'location', 'best');
    view([140, 30]);
    % title('Target and Chaser MCI Trajectories')
    if opt.saveplots
        print(fig, 'Output/Plots/Trajectory MCI.png', '-dpng', res);
    end
    
    % Chaser Acceleration in LVLH
    fig = figure('name', 'Chaser LVLH Acceleration');
    subplot(1, 3, 1)
    plot((tspan(1:size(RHOd_LVLH, 1))-t0)*TU*sec2hrs, RHOd_LVLH(:, 7)*DU/TU^2, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
    hold on
    plot((tspan-t0)*TU*sec2hrs, acc(:, 1)*DU/TU^2, 'color', '#4195e8', 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$\rho_r \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('Desired', 'Actual', 'location', 'best')
    grid on
    subplot(1, 3, 2)
    plot((tspan(1:size(RHOd_LVLH, 1))-t0)*TU*sec2hrs, RHOd_LVLH(:, 8)*DU/TU^2, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
    hold on
    plot((tspan-t0)*TU*sec2hrs, acc(:, 2)*DU/TU^2, 'color', '#4195e8', 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$\rho_\theta \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('Desired', 'Actual', 'location', 'best')
    grid on
    subplot(1, 3, 3)
    plot((tspan(1:size(RHOd_LVLH, 1))-t0)*TU*sec2hrs, RHOd_LVLH(:, 9)*DU/TU^2, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
    hold on
    plot((tspan-t0)*TU*sec2hrs, acc(:, 3)*DU/TU^2, 'color', '#4195e8', 'LineWidth', 1.5)
    xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$\rho_h \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('Desired', 'Actual', 'location', 'best')
    grid on
    if opt.saveplots
        print(fig, 'Output/Plots/Acceleration LVLH.png', '-dpng', res);
    end

    if ~sim_from_mc
        % Thrust Direction Components - xb
        fig = figure('name', "Body Thrust Direction Components");
        plot((tspan_ctrl-t0)*TU*sec2hrs, xb_AOCS_stack(:, 1), 'LineWidth', 1.5)
        hold on
        plot((tspan_ctrl-t0)*TU*sec2hrs, xb_AOCS_stack(:, 2), 'LineWidth', 1.5)
        plot((tspan_ctrl-t0)*TU*sec2hrs, xb_AOCS_stack(:, 3), 'LineWidth', 1.5)
        xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$x_{b,i}$', 'interpreter', 'latex', 'fontsize', 12)
        legend('x_{b1}', 'x_{b2}', 'x_{b3}', 'fontsize', 10, 'location', 'best')
        grid on
        if opt.saveplots
            print(fig, 'Output/Plots/Body Thrust Direction Components.png', '-dpng', res);
        end

        % Commanded Angular Velocities
        fig = figure('name', "Commanded Angular Velocities");
        plot((tspan_ctrl-t0)*TU*sec2hrs, omega_c_rt_stack(:, 1)/TU, 'LineWidth', 1.5)
        hold on
        plot((tspan_ctrl-t0)*TU*sec2hrs, omega_c_rt_stack(:, 2)/TU, 'LineWidth', 1.5)
        plot((tspan_ctrl-t0)*TU*sec2hrs, omega_c_rt_stack(:, 3)/TU, 'LineWidth', 1.5)
        xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$\omega_{c,i} \ [rad/s]$', 'interpreter', 'latex', 'fontsize', 12)
        legend('\omega_{c1}', '\omega_{c2}', '\omega_{c3}', 'fontsize', 10, 'location', 'best')
        grid on
        if opt.saveplots
            print(fig, 'Output/Plots/Commanded Angular Velocities.png', '-dpng', res);
        end
        
        % Commanded Angular Accelerations
        fig = figure('name', "Commanded Angular Accelerations");
        plot((tspan_ctrl-t0)*TU*sec2hrs, omegadot_c_rt_stack(:, 1)/TU^2, 'LineWidth', 1.5)
        hold on
        plot((tspan_ctrl-t0)*TU*sec2hrs, omegadot_c_rt_stack(:, 2)/TU^2, 'LineWidth', 1.5)
        plot((tspan_ctrl-t0)*TU*sec2hrs, omegadot_c_rt_stack(:, 3)/TU^2, 'LineWidth', 1.5)
        xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$\dot{\omega}_{c,i} \ [rad/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
        legend('\omega_{c1}', '\omega_{c2}', '\omega_{c3}', 'fontsize', 10, 'location', 'best')
        grid on
        if opt.saveplots
            print(fig, 'Output/Plots/Commanded Angular Accelerations.png', '-dpng', res);
        end

        
        % Natural vs AOCS Control Components
        fig = figure('name', "Natural vs AOCS Control Components");
        subplot(1, 3, 1)
        plot((tspan_ctrl-t0)*TU*sec2hrs, u_rt_AOCS_stack(:, 1)*1000*DU/TU^2, 'LineWidth', 1.5);
        hold on
        plot((tspan_ctrl-t0)*TU*sec2hrs, u_AOCS_stack(:, 1)*1000*DU/TU^2, 'LineWidth', 1.5);
        xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$[m/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
        legend('$u_{r, ideal}$', '$u_{r, AOCS}$', 'Location', 'best', 'Fontsize', 12, 'Interpreter','latex');
        grid on
        subplot(1, 3, 2)
        plot((tspan_ctrl-t0)*TU*sec2hrs, u_rt_AOCS_stack(:, 2)*1000*DU/TU^2, 'LineWidth', 1.5);
        hold on
        plot((tspan_ctrl-t0)*TU*sec2hrs, u_AOCS_stack(:, 2)*1000*DU/TU^2, 'LineWidth', 1.5);
        xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$[m/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
        legend('$u_{\theta, ideal}$', '$u_{\theta, AOCS}$', 'Location', 'best', 'Fontsize', 12, 'Interpreter','latex');
        grid on
        subplot(1, 3, 3)
        plot((tspan_ctrl-t0)*TU*sec2hrs, u_rt_AOCS_stack(:, 3)*1000*DU/TU^2, 'LineWidth', 1.5);
        hold on
        plot((tspan_ctrl-t0)*TU*sec2hrs, u_AOCS_stack(:, 3)*1000*DU/TU^2, 'LineWidth', 1.5);
        xlabel('$t \ [hrs]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$[m/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
        legend('$u_{h, ideal}$', '$u_{h, AOCS}$', 'Location', 'best', 'Fontsize', 12, 'Interpreter','latex');
        grid on
        if opt.saveplots
            print(fig, 'Output/Plots/Natural vs AOCS Control Components.png', '-dpng', res);
        end

    end

end


% %% Animations
% 
% show_attitude_animation = 0;
% 
% if show_attitude_animation             % show attitude evolution
%     % close all
%     clc
%     branch_animation = 1;
%     top_index = sum(indices_ctrl(1:branch_animation+1))-1;
%     fps = 10;
%     for k = 1 : M_ctrl
%         Tc = eye(4);
%         Tb = eye(4);
%         Tc(1:3, 1:3) = q2C(Q_N2C_AOCS_stack(k, 1), Q_N2C_AOCS_stack(k, 2:4)');       % rotation from MCI to Commanded
%         Tb(1:3, 1:3) = q2C(Y_ctrl(k, 14), Y_ctrl(k, 15:17)');                                  % rotation from MCI to Body
%         Tb(1:3, 1:3) = q2C(QB_LVLH(k, 1), QB_LVLH(k, 2:4)');
%         if k == 1
%             fig = figure('Name', 'Attitude Evolution');
%             % commanded = show_frame(Tc, '#349beb', 'C');
%             body = show_frame(Tb, '#fc9803', 'B');
%             N = show_frame(eye(4), '#000000', 'MCI');
%             axis([-1, 1, -1, 1, -1, 1])
%             grid on
%             fprintf('Branch: %d\n', branch_animation);
%         else
%             if k > top_index
%                 branch_animation = branch_animation + 1;
%                 top_index = sum(indices_ctrl(1:branch_animation+1))-1;
%                 fprintf('Branch: %d\n', branch_animation);
%             end
%             if rem(k, fps) == 0
%                 % update_frame(commanded, Tc);
%                 update_frame(body, Tb);
%             end
%         end
%     end
% end