function [dTCC, omega_LVLH, omegadot_LVLH, apc_LVLHt, u, rhod_LVLH,...
    rhodotd_LVLH, rhoddotd_LVLH, f, f_norm, kp_out, k_type] = HybridPredictiveControl(t, TCC, EarthPPsMCI, SunPPsMCI, muM, muE, muS, ...
    MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, tf, ppXd, DU, TU, checkTimes, Delta_t, N_inner_integration, filepath)
% Description: ...


persistent kp verifiedTimes index stopPrediction u_max

% Initialize Persistent Variables
if t == t0
    verifiedTimes = zeros(length(checkTimes),1);
    index = 1;
    stopPrediction = 0;
    u_max = 0;
end

debug = 0;


% ----- Natural Relative Motion ----- %

% Initialize State Derivative
dTCC = zeros(13, 1); 
    
% Retrieve Data from Input
MEEt = TCC(1:6)';
RHO_LVLH = TCC(7:12)';
x7 = TCC(13);

% Retrieve RHO State Variables
rho_LVLH = RHO_LVLH(1:3);
rhodot_LVLH = RHO_LVLH(4:6);

% Retrieve Target State in MCI
COEt = MEE2COE(MEEt')';
Xt_MCI = COE2rvPCI(COEt', muM)';
rt_MCI = Xt_MCI(1:3);

% Convert RHO state into MCI
RHO_MCI = rhoLVLH2MCI(RHO_LVLH, Xt_MCI, t, EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM);

% Compute Chaser State in MCI
Xc_MCI = Xt_MCI + RHO_MCI;
COEc = rvPCI2COE(Xc_MCI', muM)';
MEEc = COE2MEE(COEc')';
rc_MCI = Xc_MCI(1:3);

% Target's Third, Fourth Body and Moon Harmonics Perturbing Accelerations
a34Bt = ThirdFourthBody(MEEt, t, EarthPPsMCI, SunPPsMCI, muE, muS);
aG_Mt = MoonHarmPerts(MEEt, MoonPPsECI, t, muM, deltaE, psiM, deltaM);
apt_LVLHt = a34Bt + aG_Mt;

% Chaser's Third, Fourth Body and Moon Harmonics Perturbing Accelerations
a34Bc = ThirdFourthBody(MEEc, t, EarthPPsMCI, SunPPsMCI, muE, muS);
aG_Mc = MoonHarmPerts(MEEc, MoonPPsECI, t, muM, deltaE, psiM, deltaM);
apc_LVLHc = a34Bc + aG_Mc;

% Convert Perturbating Accelerations into MCI
[R_MCI2LVLHt, ~] = get_rotMCI2LVLH(Xt_MCI, t, EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM);
[R_LVLHc2MCI, ~] = get_rotLVLH2MCI(Xc_MCI, t, EarthPPsMCI, SunPPsMCI, MoonPPsECI, muE, muS, deltaE, psiM, deltaM);

apc_MCI = R_LVLHc2MCI*apc_LVLHc;

% Compute Angular Velocity of LVLH wrt MCI
[~, ~, incl, ~, omega, nu] = S2C(COEt');
[~, ~, incl_dot, Omega_dot, omega_dot, nu_dot] = S2C(get_COEdots(COEt', Xt_MCI', apt_LVLHt));

theta_t = omega + nu;
theta_t_dot = omega_dot + nu_dot;

omega_LVLH = [Omega_dot*sin(incl)*sin(theta_t) + incl_dot*cos(theta_t); ...
              Omega_dot*sin(incl)*cos(theta_t) - incl_dot*sin(theta_t); ...
              Omega_dot*cos(incl) + theta_t_dot];

% Compute Angular Acceleration of LVLH wrt MCI
omegadot_LVLH = ppsval(omegadotPPsLVLH, t);

% Rotate the Necessary Vectors
rt_LVLH = R_MCI2LVLHt * rt_MCI;     % this is -r_t2M^(LVLH)
apc_LVLHt = R_MCI2LVLHt * apc_MCI;

% MEE Propagation Quantities
eta = get_eta(MEEt);
G = get_G(MEEt, muM, eta);

% Relative Motion Propagation Quantities
rt = norm(rt_LVLH);
rc = norm(rc_MCI);
q = (dot(rho_LVLH, rho_LVLH) + 2*dot(rho_LVLH, rt_LVLH))/rt^2;

f = -2*cross(omega_LVLH, rhodot_LVLH) - cross(omegadot_LVLH, rho_LVLH) - cross(omega_LVLH, cross(omega_LVLH, rho_LVLH)) + muM/rt^3*((q*(2+q+(1+q)^(1/2)))/((1+q)^(3/2)*((1+q)^(1/2)+1)))*rt_LVLH - muM/rc^3*rho_LVLH + apc_LVLHt - apt_LVLHt;
f_norm = norm(f);


    
% ----- Hybrid Predictive Control ----- %

% Retrieve Reference Trajectory
RHOd_LVLH = ppsval(ppXd, t);
rhod_LVLH = RHOd_LVLH(1:3);
rhodotd_LVLH = RHOd_LVLH(4:6);
rhoddotd_LVLH = RHOd_LVLH(7:9);

zeta = 1;       % condition of critical damping

g0 = 9.80665;
u_limit =  5e-5*g0/(1000*DU)*TU^2;
threshold = 0.9*u_limit;

% ADJUSTMENT
index = min(index, length(checkTimes));

% Predictive Forward Propagation
if t > checkTimes(index) && verifiedTimes(index) == 0 && stopPrediction == 0
    
    % Define Intermediate Time Domain
    t0_int = t;
    tf_int = t + Delta_t;

    % Perform the Propagation
    [tspan, ControlledRelativeState] = odeHamHPC(@(t, state) NaturalFeedbackControl(t, state, EarthPPsMCI, SunPPsMCI, muM, ...
        muE, muS, MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, ppXd, kp, DU, TU), ...
        [t0_int, tf_int], TCC, N_inner_integration);

    % Predictive Propagation Post-Processing
    len = length(tspan);
    u_future = zeros(3,len);
    u_norms = zeros(len,1);
    for s = 1 : len
        [~, ~, ~, ~, u_future(:,s), ...
            ~, ~, ~, ~] = NaturalFeedbackControl(tspan(s), ...
            ControlledRelativeState(s,:), EarthPPsMCI, SunPPsMCI, muM, ...
            muE, muS, MoonPPsECI, deltaE, psiM, deltaM, omegadotPPsLVLH, t0, ppXd, kp, DU, TU);
        u_norms(s) = norm(u_future(:,s));
    end
    
    if debug
        
        close all
        figure('name', 'Control Debugging');

        % Visualize Control Components
        subplot(1, 2, 1)
        u1 = plot((tspan-t0)*TU/3600, u_future(1,:)*1000*DU/TU^2, 'LineWidth', 1.5);
        hold on
        u2 = plot((tspan-t0)*TU/3600, u_future(2,:)*1000*DU/TU^2, 'LineWidth', 1.5);
        u3 = plot((tspan-t0)*TU/3600, u_future(3,:)*1000*DU/TU^2, 'LineWidth', 1.5);
        ulim = plot((tspan-t0)*TU/3600, u_limit*ones(len, 1)*1000*DU/TU^2, 'r--', 'LineWidth', 1.2);
        xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$[m/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
        title('Control Components')
        grid on
        legend([u1, u2, u3, ulim], '$u_r$', '$u_{\theta}$', '$u_h$', '$u_{max}$','Location', 'best', 'Fontsize', 12, 'Interpreter','latex');
        
        % Visualize Control Norm
        subplot(1, 2, 2)
        p1 = plot((tspan-t0)*TU/3600, u_norms*1000*DU/TU^2, 'Color', '#4195e8', 'LineWidth', 1.5);
        hold on
        p2 = plot((tspan-t0)*TU/3600, u_limit*ones(len, 1)*1000*DU/TU^2, 'r--', 'LineWidth', 1.2);
        xlabel('$t \ [hours]$', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$[m/s^2]$', 'interpreter', 'latex', 'fontsize', 12)
        title('Control Norm')
        grid on
        legend([p1, p2], '$|u|$', '$u_{max}$','Location', 'best', 'Fontsize', 12, 'Interpreter', 'latex');
    
        pause
    end

    [crosses, indexes, ~] = detectCross(u_norms, threshold);
    if crosses <= 1         % this is a cross from above, that's ok
        verifiedTimes(index) = 1;
        index = index +1;
    elseif crosses >= 2     % this means we have at least a cross from below, i. e. the threshold has been overcome by the feedback reaction
        stopPrediction = 1;
        stopSaturationTime = t;
        save(filepath, "stopSaturationTime")
    end

    if isempty(indexes)
        ind = 1;
    else
        ind = indexes(1);
    end

    [u_max, ~] = max(u_norms(ind+1:end));   % first values are excluded, because it's in saturation ( therefore u_norms(1) == u_limit )

end


if u_max < threshold    % compute kp to have saturation
    a = -f + rhoddotd_LVLH;
    b = 2*(rhodot_LVLH - rhodotd_LVLH);
    c = (rho_LVLH - rhod_LVLH);

    A = norm(c)^2;
    B = 2*dot(b,c);
    C = norm(b)^2 - 2*dot(a,c);
    D = -2*dot(a,b);
    E = norm(a)^2 - u_limit^2;

    syms K              % symbolic Kp

    eqn = A*K^4 + B*K^3 + C*K^2 + D*K + E == 0;

    solutions = double(solve(eqn, K));


    % Only real positive values of kp are chosen. If imaginary part
    % is negligible, it's put to 0.
    % kp_candidates = solutions.^2;
    kp_set = [];
    for cnt = 1 : length(solutions)

        if isreal(solutions(cnt)) && solutions(cnt)>0
            kp_set = [kp_set ( solutions(cnt) )^2];
        end

    end


    % I choose the minimum positive real value of kp s.t.
    % |u| = u_lim, if it exists
    if isempty(kp_set)
        % Actual kp is imaginary, older kp is kept
        k_type = 1j;
    else
        [kp, k_type] = min(kp_set);
    end

else
    
    k_type = NaN;   % threshold overcome

end

kp_out = kp;


% Define Gain Matrices
Kp = kp*eye(3,3);
kd = 2*zeta*sqrt(kp);
Kd = kd*eye(3,3);

% Compute the Control
u = -f + rhoddotd_LVLH - Kd*(rhodot_LVLH-rhodotd_LVLH) - Kp*(rho_LVLH -rhod_LVLH);

% Compute Mass Ratio Derivative
c = 30/DU*TU;           % effective exhaust velocity = 30 km/s
x7_dot = - x7*norm(u)/c;



% Assign State Derivatives
dTCC(1:6) = G*apt_LVLHt;
dTCC(6) = dTCC(6) + sqrt(muM/MEEt(1)^3)*eta^2;
dTCC(7:9) = rhodot_LVLH;
dTCC(10:12) = f + u;
dTCC(13) = x7_dot;


% Clock for the Integration
global pbar opt

if opt.show_progress
    Day = 86400;  % seconds in a day
    Hour = 3600;  % seconds in an hour
    Min = 60;     % seconds in a minute
    tDAY = floor((t - t0) * TU / Day);      % calculate the elapsed time components
    tHR = floor(((t - t0) * TU - tDAY * Day) / Hour);
    tMIN = floor(((t - t0) * TU - tDAY * Day - tHR * Hour) / Min);
    timeStr = sprintf('Time Elapsed: %02d days, %02d hrs, %02d mins', tDAY, tHR, tMIN);     % create a string for the time
    waitbarMessage = sprintf('Progress: %.2f%%\n%s', (t-t0)/(tf-t0)*100, timeStr);      % create the waitbar message including the time and progress percentage
    waitbar((t-t0)/(tf-t0), pbar, waitbarMessage);      % update the waitbar
end

end

