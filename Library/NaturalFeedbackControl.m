function [dTCC, omega_LVLH, omegadot_LVLH, apc_LVLHt, u, rhod_LVLH,...
    rhodotd_LVLH, rhoddotd_LVLH, f_norm] = NaturalFeedbackControl(t, ...
    TCC, EarthPPsMCI, SunPPsMCI, muM, muE, muS, MoonPPsECI, deltaE, ...
    psiM, deltaM, omegadotPPsLVLH, t0, tf, ppXd, kp, u_lim, DU, TU, misalignment, clock, is_col, emergency_manoeuvre_flag, rho_0, rho_f)

persistent emergency_hysteresis

% ----- Natural Relative Motion ----- %

global pbar

% Initialize State Derivative
dTCC = zeros(13, 1); 
    
% Retrieve Data from Input

if is_col
    MEEt = TCC(1:6);
    RHO_LVLH = TCC(7:12);
else
    MEEt = TCC(1:6)';
    RHO_LVLH = TCC(7:12)';
end
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


    
% ----- Natural Feedback Control ----- %

% Retrieve Reference Trajectory
RHOd_LVLH = ppsval(ppXd, t);
rhod_LVLH = RHOd_LVLH(1:3);
rhodotd_LVLH = RHOd_LVLH(4:6);
rhoddotd_LVLH = RHOd_LVLH(7:9);

% Define Gain Matrices
Kp = kp*eye(3,3);
zeta = 1;       % condition of critical damping
kd = 2*zeta*sqrt(kp);
Kd = kd*eye(3,3);

% Compute the Nominal Control Thrust
un = -f + rhoddotd_LVLH - Kd*(rhodot_LVLH-rhodotd_LVLH) - Kp*(rho_LVLH -rhod_LVLH);
un_hat = un / norm(un);
un_norm = norm(un);
if un_norm > u_lim
    un_norm = u_lim;
end


% Apply Thrust Misalignment
alphan = atan2(un_hat(2), un_hat(1));
deltan = asin(un_hat(3));
beta = misalignment.beta;
gamma = misalignment.gamma;

ur = un_norm * (sin(gamma) * cos(beta) * sin(alphan) + sin(gamma) * sin(beta) * cos(deltan) * cos(alphan) + cos(gamma) * cos(deltan) * cos(alphan));
ut = un_norm * (-sin(gamma) * cos(beta) * cos(alphan) + sin(gamma) * sin(beta) * sin(deltan) * sin(alphan) + cos(gamma) * cos(deltan) * sin(alphan));
uh = un_norm * (-sin(gamma) * sin(beta) * cos(deltan) + cos(gamma) * sin(deltan));
u = [ur; ut; uh];

% Apply Emergency Manoeuvre
if emergency_manoeuvre_flag

    emergency_min_distance = 14.8*1e-3/DU;           % 9.8 m
    emergency_max_distance = 15.2*1e-3/DU;          % 10.2 m
    
    if isempty(emergency_hysteresis)
        if norm(rho_LVLH) < emergency_min_distance
            emergency_hysteresis = 1;
        else
            emergency_hysteresis = 0;
        end
    end
    
    if emergency_hysteresis && norm(rho_LVLH) > emergency_max_distance  
        emergency_hysteresis = 0;
        fprintf('Emergency Manoeuvre: 1 -> 0\n');
    end
    
    if ~emergency_hysteresis && norm(rho_LVLH) < emergency_min_distance
        emergency_hysteresis = 1;
        fprintf('Emergency Manoeuvre: 0 -> 1\n');
    end
    
    if norm(rho_LVLH) < emergency_min_distance
        l_hat = cross(rho_0, rho_f) / norm(cross(rho_0, rho_f));
        lambda0_hat = cross(l_hat, rho_0)/norm(rho_0);
        u = u_lim * lambda0_hat;
        fprintf('Applied Emergency Manoeuvre: dist = %.3f m\n', norm(rho_LVLH)*DU*1e3);
    elseif norm(rho_LVLH) >= emergency_min_distance && norm(rho_LVLH) <= emergency_max_distance && emergency_hysteresis
        l_hat = cross(rho_0, rho_f) / norm(cross(rho_0, rho_f));
        lambda0_hat = cross(l_hat, rho_0)/norm(rho_0);
        u = u_lim * lambda0_hat;
        fprintf('Applied Emergency Manoeuvre: dist = %.3f m\n', norm(rho_LVLH)*DU*1e3);
    end

end

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
if clock
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

