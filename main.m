%% Moon Harbour Project - Leonardo Russo

close all
clear
clc
                                                                                
addpath('Library/')
addpath('Data/')
addpath('Data/Planets/')
addpath('Data/Materials/')
addpath('Data/Ephemeris/')

sampling_time = 60;                     % seconds
include_actuation = false;
verbose = true;
misalignment_type = "oscillating";
state_perturbation_flag = true;
engine_failure_flag = true;
workspace_path = "Data/Utils/main.mat";

parfmain(sampling_time, include_actuation, verbose, misalignment_type, state_perturbation_flag, engine_failure_flag, workspace_path);


%% Visualize the Results

close all

load('Data/Utils/main.mat');

opt.saveplots = false;
opt.additional_plots = false;

% Show Numerical Results
fprintf('\nTotal Runtime: %.1f s.\n\n', runtime)
VarNames = {'r (m)', 'theta (m)', 'h (m)', 'v_r (m/s)', 'v_theta (m/s)', 'v_h (m/s)'};
DesRhoState = array2table(desState, 'VariableNames', VarNames, 'RowNames', {'Desired State'});
FinalRhoState = array2table(finalState, 'VariableNames', VarNames, 'RowNames', {'Final State'});
DeltaRhoState = array2table(deltaState, 'VariableNames', VarNames, 'RowNames', {'Delta State'});
ResultsTable = [DesRhoState; FinalRhoState; DeltaRhoState];      % concatenate the two tables vertically
disp(ResultsTable);


% Terminal Chaser State in LVLH
figure('name', 'Terminal Chaser Trajectory in LVLH Space')
C_LVLH_T = DrawTrajLVLH3D(RHO_LVLH(M_ctrl_DA:end, 1:3)*DU);
Cd_LVLH_T = DrawTrajLVLH3D(RHOd_LVLH(M_ctrl_DA:end, 1:3)*DU, '#6efad2', '-.');
% vp_T = plot3(viapoints_T(:, 1)*DU, viapoints_T(:, 2)*DU, viapoints_T(:, 3)*DU, 'color', 'r', 'linestyle', 'none', 'marker', '.', 'markersize', 15);
title('Terminal Chaser LVLH Trajectory')
% legend([C_LVLH_T, Cd_LVLH_T, vp_T], {'Chaser Trajectory', 'Reference Trajectory', 'Terminal Via Points'}, 'location', 'best')
legend([C_LVLH_T, Cd_LVLH_T], {'Chaser Trajectory', 'Reference Trajectory'}, 'location', 'best')
if opt.saveplots
    saveas(gcf, strcat('Output/Plots/Trajectory Terminal LVLH.jpg'))
end


% Chaser LVLH State Components
figure('name', 'Chaser LVLH State Components')
subplot(2, 3, 1)
plot((tspan(1:size(RHOd_LVLH, 1))-t0)*TU*sec2hrs, RHOd_LVLH(:, 1)*DU, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
hold on
plot((tspan-t0)*TU*sec2hrs, RHO_LVLH(:, 1)*DU, 'color', '#4195e8', 'LineWidth', 1.5)
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\rho_r \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('Desired', 'Actual', 'fontsize', 10, 'location', 'best')
grid on
subplot(2, 3, 2)
plot((tspan(1:size(RHOd_LVLH, 1))-t0)*TU*sec2hrs, RHOd_LVLH(:, 2)*DU, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
hold on
plot((tspan-t0)*TU*sec2hrs, RHO_LVLH(:, 2)*DU, 'color', '#4195e8', 'LineWidth', 1.5)
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\rho_\theta \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('Desired', 'Actual', 'fontsize', 10, 'location', 'best')
grid on
subplot(2, 3, 3)
plot((tspan(1:size(RHOd_LVLH, 1))-t0)*TU*sec2hrs, RHOd_LVLH(:, 3)*DU, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
hold on
plot((tspan-t0)*TU*sec2hrs, RHO_LVLH(:, 3)*DU, 'color', '#4195e8', 'LineWidth', 1.5)
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\rho_h \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('Desired', 'Actual', 'fontsize', 10, 'location', 'best')
grid on
subplot(2, 3, 4)
plot((tspan(1:size(RHOd_LVLH, 1))-t0)*TU*sec2hrs, RHOd_LVLH(:, 4)*DU/TU, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
hold on
plot((tspan-t0)*TU*sec2hrs, RHO_LVLH(:, 4)*DU/TU, 'color', '#4195e8', 'LineWidth', 1.5)
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\dot{\rho}_r \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('Desired', 'Actual', 'fontsize', 10, 'location', 'best')
grid on
subplot(2, 3, 5)
plot((tspan(1:size(RHOd_LVLH, 1))-t0)*TU*sec2hrs, RHOd_LVLH(:, 5)*DU/TU, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
hold on
plot((tspan-t0)*TU*sec2hrs, RHO_LVLH(:, 5)*DU/TU, 'color', '#4195e8', 'LineWidth', 1.5)
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\dot{\rho}_\theta \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('Desired', 'Actual', 'fontsize', 10, 'location', 'best')
grid on
subplot(2, 3, 6)
plot((tspan(1:size(RHOd_LVLH, 1))-t0)*TU*sec2hrs, RHOd_LVLH(:, 6)*DU/TU, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
hold on
plot((tspan-t0)*TU*sec2hrs, RHO_LVLH(:, 6)*DU/TU, 'color', '#4195e8', 'LineWidth', 1.5)
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\dot{\rho}_h \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
legend('Desired', 'Actual', 'fontsize', 10, 'location', 'best')
grid on
if opt.saveplots
    saveas(gcf, strcat('Output/Plots/State LVLH Components.jpg'))
end


% Distance and Velocity Norms
figure('Name', 'Distance and Velocity Norms')
subplot(1, 2, 1)
plot((tspan-t0)*TU*sec2hrs, dist*DU, 'color', '#4195e8', 'LineWidth', 1.5)
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$|\rho| \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
title('Relative Distance')
grid on
subplot(1, 2, 2)
plot((tspan-t0)*TU*sec2hrs, vel*DU/TU, 'color', '#4195e8', 'LineWidth', 1.5)
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$|v| \ [km/s]$', 'interpreter', 'latex', 'fontsize', 12)
title('Relative Velocity')
grid on
if opt.saveplots
    saveas(gcf, strcat('Output/Plots/Relative Distance and Velocity.jpg'))
end


% Control Norm
figure('name', 'Control Thrust')
p1 = plot((tspan-t0)*TU*sec2hrs, u_norms*1000*DU/TU^2, 'Color', '#4195e8', 'LineWidth', 1.5);
hold on
p2 = plot((tspan-t0)*TU*sec2hrs, u_limit*ones(length(tspan), 1), 'r--', 'LineWidth', 1.2);
p3 = plot((tspan-t0)*TU*sec2hrs, f_norms*1000*DU/TU^2, 'Color', '#93faad', 'LineWidth', 1.5);
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$[m/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
title('Control Norm')
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
xlim([(t0_T-t0)*TU*sec2hrs, (tf-t0)*TU*sec2hrs]);   % set the x and y limits for the zoomed plot based on the final part of the data
% ylim([-u_limit(1)/2, u_limit(1)/2]);
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 10)
ylabel('$[m/s^2]$', 'interpreter', 'latex', 'fontsize', 10)
if opt.saveplots
    saveas(gcf, strcat('Output/Plots/Control Norm.jpg'))
end



% Control Components
figure('name', 'Control Thrust Components')
u1 = plot((tspan-t0)*TU*sec2hrs, u(:, 1)*1000*DU/TU^2, 'LineWidth', 1.5);
hold on
u2 = plot((tspan-t0)*TU*sec2hrs, u(:, 2)*1000*DU/TU^2, 'LineWidth', 1.5);
u3 = plot((tspan-t0)*TU*sec2hrs, u(:, 3)*1000*DU/TU^2, 'LineWidth', 1.5);
ulim = plot((tspan-t0)*TU*sec2hrs, u_limit*ones(length(tspan), 1), 'r--', 'LineWidth', 1.2);
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$[m/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
title('Control Components')
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
xlim([(t0_T-t0)*TU*sec2hrs, (tf-t0)*TU*sec2hrs]); % set the x and y limits to focus on the final part of the plot
% ylim([-u_limit(1)/2, u_limit(1)/2]);
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 10)
ylabel('$[m/s^2]$', 'interpreter', 'latex', 'fontsize', 10)
if opt.saveplots
    saveas(gcf, strcat('Output/Plots/Control Components.jpg'))
end


% Body and Commanded Attitude Quaternions
figure('name', "Body and Commanded Attitude")
subplot(1, 2, 1)
plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 14), 'LineWidth', 1.5)
hold on
plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 15), 'LineWidth', 1.5)
plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 16), 'LineWidth', 1.5)
plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 17), 'LineWidth', 1.5)
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$q_i$', 'interpreter', 'latex', 'fontsize', 12)
legend('q_{b0}', 'q_{b1}', 'q_{b2}', 'q_{b3}', 'fontsize', 10, 'location', 'best')
grid on
subplot(1, 2, 2)
plot((tspan_ctrl-t0)*TU*sec2hrs, Q_N2C_AOCS_stack(:, 1), 'LineWidth', 1.5)
hold on
plot((tspan_ctrl-t0)*TU*sec2hrs, Q_N2C_AOCS_stack(:, 2), 'LineWidth', 1.5)
plot((tspan_ctrl-t0)*TU*sec2hrs, Q_N2C_AOCS_stack(:, 3), 'LineWidth', 1.5)
plot((tspan_ctrl-t0)*TU*sec2hrs, Q_N2C_AOCS_stack(:, 4), 'LineWidth', 1.5)
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$q_i$', 'interpreter', 'latex', 'fontsize', 12)
legend('q_{c0}', 'q_{c1}', 'q_{c2}', 'q_{c3}', 'fontsize', 10, 'location', 'best')
grid on  
if opt.saveplots
    saveas(gcf, strcat('Output/Plots/Body and Commanded Quaternions.jpg'))
end


% Error Quaternions
figure('name', "Error Quaternions")
plot((tspan_ctrl-t0)*TU*sec2hrs, qe0, 'LineWidth', 1.5)
hold on
plot((tspan_ctrl-t0)*TU*sec2hrs, qe(:, 1), 'LineWidth', 1.5)
plot((tspan_ctrl-t0)*TU*sec2hrs, qe(:, 2), 'LineWidth', 1.5)
plot((tspan_ctrl-t0)*TU*sec2hrs, qe(:, 3), 'LineWidth', 1.5)
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$q_i$', 'interpreter', 'latex', 'fontsize', 12)
legend('q_{e0}', 'q_{e1}', 'q_{e2}', 'q_{e3}', 'fontsize', 10, 'location', 'best')
grid on
if opt.saveplots
    saveas(gcf, strcat('Output/Plots/Error Quaternions.jpg'))
end


% Body and Wheels Angular Velocities
figure('name', "Body and Wheels Angular Velocities")
subplot(1, 2, 1)
plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 18)/TU, 'LineWidth', 1.5)
hold on
plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 19)/TU, 'LineWidth', 1.5)
plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 20)/TU, 'LineWidth', 1.5)
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\omega_i \, [rad/s]$', 'interpreter', 'latex', 'fontsize', 12)
legend('\omega_{1}', '\omega_{2}', '\omega_{3}', 'fontsize', 10, 'location', 'best')
grid on
subplot(1, 2, 2)
plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 21)/TU, 'LineWidth', 1.5)
hold on
plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 22)/TU, 'LineWidth', 1.5)
plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 23)/TU, 'LineWidth', 1.5)
plot((tspan_ctrl-t0)*TU*sec2hrs, Y_ctrl(:, 24)/TU, 'LineWidth', 1.5)
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\omega_{si} \, [rad/s]$', 'interpreter', 'latex', 'fontsize', 12)
legend('\omega_{s1}', '\omega_{s2}', '\omega_{s3}', '\omega_{s4}', 'fontsize', 10, 'location', 'best')
grid on
if opt.saveplots
    saveas(gcf, strcat('Output/Plots/Body and Wheels Angular Velocities.jpg'))
end


% Commanded Torque
figure('name', 'Commanded Torque')
plot((tspan_ctrl-t0)*TU*sec2hrs, Tc(:, 1)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
hold on
plot((tspan_ctrl-t0)*TU*sec2hrs, Tc(:, 2)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
plot((tspan_ctrl-t0)*TU*sec2hrs, Tc(:, 3)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$T_{c,i} \ [Nm]$', 'interpreter', 'latex', 'fontsize', 12)
legend('T_{c1}', 'T_{c2}', 'T_{c3}', 'fontsize', 10, 'location', 'best')
grid on
if opt.saveplots
    saveas(gcf, strcat('Output/Plots/Commanded Torque.jpg'))
end

if opt.include_actuation
    % Actual Torque
    figure('name', 'Actual Torque')
    plot((tspan_ctrl-t0)*TU*sec2hrs, Ta(:, 1)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
    hold on
    plot((tspan_ctrl-t0)*TU*sec2hrs, Ta(:, 2)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
    plot((tspan_ctrl-t0)*TU*sec2hrs, Ta(:, 3)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
    xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$T_{a,i} \ [Nm]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('T_{a1}', 'T_{a2}', 'T_{a3}', 'fontsize', 10, 'location', 'best')
    grid on
    if opt.saveplots
        saveas(gcf, strcat('Output/Plots/Actual Torque.jpg'))
    end
    
    
    % Actual vs Commanded Torque
    figure('name', 'Actual vs Commanded Torque')
    subplot(1, 3, 1)
    plot((tspan_ctrl-t0)*TU*sec2hrs, Ta(:, 1)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
    hold on
    plot((tspan_ctrl-t0)*TU*sec2hrs, Tc(:, 1)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
    xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$T_{1} \ [Nm]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('T_{a1}', 'T_{c1}', 'fontsize', 10, 'location', 'best')
    grid on
    subplot(1, 3, 2)
    plot((tspan_ctrl-t0)*TU*sec2hrs, Ta(:, 2)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
    hold on
    plot((tspan_ctrl-t0)*TU*sec2hrs, Tc(:, 2)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
    xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$T_{2} \ [Nm]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('T_{a2}', 'T_{c2}', 'fontsize', 10, 'location', 'best')
    grid on
    subplot(1, 3, 3)
    plot((tspan_ctrl-t0)*TU*sec2hrs, Ta(:, 3)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
    hold on
    plot((tspan_ctrl-t0)*TU*sec2hrs, Tc(:, 3)*1e6*DU^2/TU^2*MU, 'LineWidth', 1.5)
    xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$T_{3} \ [Nm]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('T_{a3}', 'T_{c3}', 'fontsize', 10, 'location', 'best')
    grid on
    if opt.saveplots
        saveas(gcf, strcat('Output/Plots/Actual vs Commanded Torque.jpg'))
    end
end


% Gammas and Betas
if misalignment_type ~= "null"
    figure('name', "Misalignment Angles")
    subplot(1, 2, 1)
    plot((tspan_ctrl-t0)*TU*sec2hrs, rad2deg(betas), 'LineWidth', 1.5)
    xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$\beta \ [deg]$', 'interpreter', 'latex', 'fontsize', 12)
    grid on
    subplot(1, 2, 2)
    plot((tspan_ctrl-t0)*TU*sec2hrs, rad2deg(gammas), 'LineWidth', 1.5)
    xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$\gamma \ [deg]$', 'interpreter', 'latex', 'fontsize', 12)
    grid on
    if opt.saveplots
        saveas(gcf, strcat('Output/Plots/Beta and Gamma angles.jpg'))
    end
end


if opt.additional_plots

    % Chaser State in LVLH
    figure('name', 'Chaser Trajectory in LVLH Space')
    C_LVLH = DrawTrajLVLH3D(RHO_LVLH(:, 1:3)*DU);
    Cd_LVLH = DrawTrajLVLH3D(RHOd_LVLH(:, 1:3)*DU, '#6efad2', '-.');
    % vp_T = plot3(viapoints_T(:, 1)*DU, viapoints_T(:, 2)*DU, viapoints_T(:, 3)*DU, 'color', 'r', 'linestyle', 'none', 'marker', '.', 'markersize', 15);
    vp_DA = plot3(viapoints_DA(:, 1)*DU, viapoints_DA(:, 2)*DU, viapoints_DA(:, 3)*DU, 'color', 'b', 'linestyle', 'none', 'marker', '.', 'markersize', 15);
    title('Chaser LVLH Trajectory')
    % legend([C_LVLH, Cd_LVLH, vp_T, vp_DA], {'Chaser Trajectory', 'Reference Trajectory', 'Terminal Via Points', 'Direct Approach Via Points'}, 'location', 'best')
    legend([C_LVLH, Cd_LVLH, vp_DA], {'Chaser Trajectory', 'Reference Trajectory', 'Direct Approach Via Points'}, 'location', 'best')
    if opt.saveplots
        saveas(gcf, strcat('Output/Trajectory LVLH.jpg'))
    end

    % Target, Chaser and Reference Chaser Trajectories in MCI
    figure('name', 'Trajectory in MCI Space')
    T = DrawTrajMCI3D(Xt_MCI(:, 1:3)*DU, '#d1d1d1', '-.');
    C = DrawTrajMCI3D(Xc_MCI(:, 1:3)*DU);
    legend([T, C], {'Target Trajectory', 'Chaser Trajectory'}, 'location', 'best');
    view([140, 30]);
    title('Target and Chaser MCI Trajectories')
    if opt.saveplots
        saveas(gcf, strcat('Output/Trajectory MCI.jpg'))
    end
    
    % Chaser Acceleration in LVLH
    figure('name', 'Chaser LVLH Acceleration')
    subplot(1, 3, 1)
    plot((tspan(1:size(RHOd_LVLH, 1))-t0)*TU*sec2hrs, RHOd_LVLH(:, 7)*DU/TU^2, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
    hold on
    plot((tspan-t0)*TU*sec2hrs, acc(:, 1)*DU/TU^2, 'color', '#4195e8', 'LineWidth', 1.5)
    xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$\rho_r \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('Desired', 'Actual', 'fontsize', 10, 'location', 'best')
    grid on
    
    subplot(1, 3, 2)
    plot((tspan(1:size(RHOd_LVLH, 1))-t0)*TU*sec2hrs, RHOd_LVLH(:, 8)*DU/TU^2, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
    hold on
    plot((tspan-t0)*TU*sec2hrs, acc(:, 2)*DU/TU^2, 'color', '#4195e8', 'LineWidth', 1.5)
    xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$\rho_\theta \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('Desired', 'Actual', 'fontsize', 10, 'location', 'best')
    grid on
    
    subplot(1, 3, 3)
    plot((tspan(1:size(RHOd_LVLH, 1))-t0)*TU*sec2hrs, RHOd_LVLH(:, 9)*DU/TU^2, 'color', '#6efad2', 'LineStyle', '-.', 'LineWidth', 1.2)
    hold on
    plot((tspan-t0)*TU*sec2hrs, acc(:, 3)*DU/TU^2, 'color', '#4195e8', 'LineWidth', 1.5)
    xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$\rho_h \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('Desired', 'Actual', 'fontsize', 10, 'location', 'best')
    grid on

    % Thrust Direction Components - xb
    figure('name', "Body Thrust Direction Components")
    plot((tspan_ctrl-t0)*TU*sec2hrs, xb_AOCS_stack(:, 1), 'LineWidth', 1.5)
    hold on
    plot((tspan_ctrl-t0)*TU*sec2hrs, xb_AOCS_stack(:, 2), 'LineWidth', 1.5)
    plot((tspan_ctrl-t0)*TU*sec2hrs, xb_AOCS_stack(:, 3), 'LineWidth', 1.5)
    xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$x_{b,i}$', 'interpreter', 'latex', 'fontsize', 12)
    legend('x_{b1}', 'x_{b2}', 'x_{b3}', 'fontsize', 10, 'location', 'best')
    grid on
    if opt.saveplots
        saveas(gcf, strcat('Output/Plots/Body Thrust Direction Components.jpg'))
    end

    % Commanded Angular Velocities
    figure('name', "Commanded Angular Velocities")
    plot((tspan_ctrl-t0)*TU*sec2hrs, omega_c_rt_stack(:, 1)/TU, 'LineWidth', 1.5)
    hold on
    plot((tspan_ctrl-t0)*TU*sec2hrs, omega_c_rt_stack(:, 2)/TU, 'LineWidth', 1.5)
    plot((tspan_ctrl-t0)*TU*sec2hrs, omega_c_rt_stack(:, 3)/TU, 'LineWidth', 1.5)
    xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$\omega_{c,i} \ [rad/s]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('\omega_{c1}', '\omega_{c2}', '\omega_{c3}', 'fontsize', 10, 'location', 'best')
    grid on
    if opt.saveplots
        saveas(gcf, strcat('Output/Plots/Commanded Angular Velocities.jpg'))
    end
    
    
    % Commanded Angular Accelerations
    figure('name', "Commanded Angular Accelerations")
    plot((tspan_ctrl-t0)*TU*sec2hrs, omegadot_c_rt_stack(:, 1)/TU^2, 'LineWidth', 1.5)
    hold on
    plot((tspan_ctrl-t0)*TU*sec2hrs, omegadot_c_rt_stack(:, 2)/TU^2, 'LineWidth', 1.5)
    plot((tspan_ctrl-t0)*TU*sec2hrs, omegadot_c_rt_stack(:, 3)/TU^2, 'LineWidth', 1.5)
    xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$\dot{\omega}_{c,i} \ [rad/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('\omega_{c1}', '\omega_{c2}', '\omega_{c3}', 'fontsize', 10, 'location', 'best')
    grid on
    if opt.saveplots
        saveas(gcf, strcat('Output/Plots/Commanded Angular Accelerations.jpg'))
    end

    % Natural vs AOCS Control Components
    figure('name', "Natural vs AOCS Control Components")
    subplot(1, 3, 1)
    plot((tspan_ctrl-t0)*TU*sec2hrs, u_rt_AOCS_stack(:, 1)*1000*DU/TU^2, 'LineWidth', 1.5);
    hold on
    plot((tspan_ctrl-t0)*TU*sec2hrs, u_AOCS_stack(:, 1)*1000*DU/TU^2, 'LineWidth', 1.5);
    xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$[m/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('$u_{r, ideal}$', '$u_{r, AOCS}$', 'Location', 'best', 'Fontsize', 12, 'Interpreter','latex');
    grid on
    subplot(1, 3, 2)
    plot((tspan_ctrl-t0)*TU*sec2hrs, u_rt_AOCS_stack(:, 2)*1000*DU/TU^2, 'LineWidth', 1.5);
    hold on
    plot((tspan_ctrl-t0)*TU*sec2hrs, u_AOCS_stack(:, 2)*1000*DU/TU^2, 'LineWidth', 1.5);
    xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$[m/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('$u_{\theta, ideal}$', '$u_{\theta, AOCS}$', 'Location', 'best', 'Fontsize', 12, 'Interpreter','latex');
    grid on
    subplot(1, 3, 3)
    plot((tspan_ctrl-t0)*TU*sec2hrs, u_rt_AOCS_stack(:, 3)*1000*DU/TU^2, 'LineWidth', 1.5);
    hold on
    plot((tspan_ctrl-t0)*TU*sec2hrs, u_AOCS_stack(:, 3)*1000*DU/TU^2, 'LineWidth', 1.5);
    xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
    ylabel('$[m/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
    legend('$u_{h, ideal}$', '$u_{h, AOCS}$', 'Location', 'best', 'Fontsize', 12, 'Interpreter','latex');
    grid on
    if opt.saveplots
        saveas(gcf, strcat('Output/Plots/Natural vs AOCS Control Components.jpg'))
    end

end


%% Animations

show_attitude_animation = 0;

if show_attitude_animation             % show attitude evolution
    % close all
    clc
    branch_animation = 1;
    top_index = sum(indices_ctrl(1:branch_animation+1))-1;
    fps = 10;
    for k = 1 : M_ctrl
        Tc = eye(4);
        Tb = eye(4);
        Tc(1:3, 1:3) = q2C(Q_N2C_AOCS_stack(k, 1), Q_N2C_AOCS_stack(k, 2:4)');       % rotation from MCI to Commanded
        Tb(1:3, 1:3) = q2C(Y_ctrl(k, 14), Y_ctrl(k, 15:17)');                                  % rotation from MCI to Body
        Tb(1:3, 1:3) = q2C(QB_LVLH(k, 1), QB_LVLH(k, 2:4)');
        if k == 1
            figure('Name', 'Attitude Evolution');
            % commanded = show_frame(Tc, '#349beb', 'C');
            body = show_frame(Tb, '#fc9803', 'B');
            N = show_frame(eye(4), '#000000', 'MCI');
            axis([-1, 1, -1, 1, -1, 1])
            grid on
            fprintf('Branch: %d\n', branch_animation);
        else
            if k > top_index
                branch_animation = branch_animation + 1;
                top_index = sum(indices_ctrl(1:branch_animation+1))-1;
                fprintf('Branch: %d\n', branch_animation);
            end
            if rem(k, fps) == 0
                % update_frame(commanded, Tc);
                update_frame(body, Tb);
            end
        end
    end
end