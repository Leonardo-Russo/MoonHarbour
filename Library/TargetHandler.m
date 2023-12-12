function [omegaPPsLVLH, omegadotPPsLVLH] = TargetHandler(Xt_MCI, COEt, MEEt, tspan, ...
          EarthPPsMCI, SunPPsMCI, MoonPPsECI, deltaE, psiM, deltaM, muE, muS)


% Recall Global Variables
global DU TU muM Rm

Day = 86400;

% Initialize Local Variables
omega_LVLH = zeros(length(tspan), 3);

% Compute Local Variables
pbar = waitbar(0, 'Performing the Target Trajectory Interpolation');

for i = 1 : length(tspan)

    ti = tspan(i);      % get current epoch

    % Compute Third, Fourth Body and Moon Harmonics Perturbing Accelerations
    a34B = ThirdFourthBody(MEEt(i, :)', ti, EarthPPsMCI, SunPPsMCI, muE, muS);
    aG_M = MoonHarmPerts(MEEt(i, :)', MoonPPsECI, ti, muM, deltaE, psiM, deltaM);
    ap_LVLHt = a34B + aG_M;

    % Retrieve COE Values
    [~, ~, incl, ~, omega, nu] = S2C(COEt(i, :));
    theta_t = omega + nu;

    % Compute COE Derivatives
    COEtdots = get_COEdots(COEt(i, :), Xt_MCI(i, :), ap_LVLHt);
    [~, ~, incl_dot, Omega_dot, omega_dot, nu_dot] = S2C(COEtdots);
    theta_t_dot = omega_dot + nu_dot;

    % Compute LVLH Angular Velocity wrt MCI
    omega_LVLH(i, :) = [Omega_dot*sin(incl)*sin(theta_t) + incl_dot*cos(theta_t); ...
                        Omega_dot*sin(incl)*cos(theta_t) - incl_dot*sin(theta_t); ...
                        Omega_dot*cos(incl) + theta_t_dot]';

    waitbarMessage = sprintf('Target Interpolation Progress: %.2f%%\n', i/length(tspan)*100);
    waitbar(i/length(tspan), pbar, waitbarMessage);      % update the waitbar

end

% % Show the Evolution of omegaLVLH
% % Visualize the Evolution of omega_r_LVLH
% figure('name', 'Evolution of omega_LVLH_r')
% plot((tspan - tspan(1))*TU/Day, omega_LVLH(:, 1)/TU)
% title('Evolution of omega^{(LVLH)}_r')
% xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 12)
% ylabel('$\omega_r^{(LVLH)} \ [rad/s]$', 'Interpreter','latex', 'FontSize', 12)
% 
% % Visualize the Evolution of omega_t_LVLH
% figure('name', 'Evolution of omega_LVLH_t')
% plot((tspan - tspan(1))*TU/Day, omega_LVLH(:, 2)/TU)
% title('Evolution of omega^{(LVLH)}_t')
% xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 12)
% ylabel('$\omega_{\theta}^{(LVLH)} \ [rad/s]$', 'Interpreter','latex', 'FontSize', 12)
% 
% % Visualize the Evolution of omega_h_LVLH
% figure('name', 'Evolution of omega_LVLH_h')
% plot((tspan - tspan(1))*TU/Day, omega_LVLH(:, 3)/TU)
% title('Evolution of omega^{(LVLH)}_h')
% xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 12)
% ylabel('$\omega_h^{(LVLH)} \ [rad/s]$', 'Interpreter','latex', 'FontSize', 12)

% Perform the Interpolations
omegaPPsLVLH = get_statePP(tspan, omega_LVLH);

omegadot_rPPs = fnder(omegaPPsLVLH(1), 1);
omegadot_tPPs = fnder(omegaPPsLVLH(2), 1);
omegadot_hPPs = fnder(omegaPPsLVLH(3), 1);
omegadotPPsLVLH = [omegadot_rPPs; omegadot_tPPs; omegadot_hPPs];

close(pbar);

end