function [dY, omega_LVLH, omegadot_LVLH, apc_LVLHt, u, rhod_LVLH,...
    rhodotd_LVLH, rhoddotd_LVLH, f_norm] = AOCS(t, Y, EarthPPsMCI, SunPPsMCI, muM, muE, muS, MoonPPsECI, deltaE, ...
    psiM, deltaM, omegadotPPsLVLH, t0, tf, ppXd, kp, DU, TU, omega_cPPs, omegadot_cPPs, R_N2C_PPs, sign_qe0_0, misalignment, clock)

% -------------------- Orbital Control -------------------- %

% ----- Natural Relative Motion ----- %

global pbar

% Initialize State Derivative
dY = zeros(24, 1); 
    
% Retrieve Data from Input
MEEt = Y(1:6)';
RHO_LVLH = Y(7:12)';
x7 = Y(13);

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

% Apply Thrust Misalignment
alphan = atan2(un_hat(2), un_hat(1));
deltan = asin(un_hat(3));
beta = misalignment.beta;
gamma = misalignment.gamma;

ur = un_norm * (sin(gamma) * cos(beta) * sin(alphan) + sin(gamma) * sin(beta) * cos(deltan) * cos(alphan) + cos(gamma) * cos(deltan) * cos(alphan));
ut = un_norm * (-sin(gamma) * cos(beta) * cos(alphan) + sin(gamma) * sin(beta) * sin(deltan) * sin(alphan) + cos(gamma) * cos(deltan) * sin(alphan));
uh = un_norm * (-sin(gamma) * sin(beta) * cos(deltan) + cos(gamma) * sin(deltan));
u = [ur; ut; uh];

% Compute Mass Ratio Derivative
c = 30/DU*TU;           % effective exhaust velocity = 30 km/s
x7_dot = - x7*norm(u)/c;



% Assign State Derivatives
dY(1:6) = G*apt_LVLHt;
dY(6) = dY(6) + sqrt(muM/MEEt(1)^3)*eta^2;
dY(7:9) = rhodot_LVLH;
dY(10:12) = f + u;
dY(13) = x7_dot;



% -------------------- Attitude Control -------------------- %

% Retrieve Attitude State Variables
xb = Y(14:17)';
w = Y(18:20)';
omegas = Y(21:24)';

qb0 = xb(1);
qb = xb(2:4);

% Define Chaser Parameters
Jc = [900, 50, -100;...
      50, 1100, 150;...
      -100, 150, 1250];     % kg m^2

% Gain Parameters
omega_n = 0.03;     % rad/s
omega_n = 10;
xi = 1;
c1 = 2 * omega_n^2;
c2 = 2*xi*omega_n/c1;
invA = c1 * eye(3);
B = c2 * eye(3);

% Retrieve Commanded Attitude
wc = ppsval(omega_cPPs, t);
wc_dot = ppsval(omegadot_cPPs, t);
wd = w - wc;                            % compute desired attitude

% Compute Error Quaternions
R_N2C = rotppsval(R_N2C_PPs, t);
[qc0, qc] = C2q(R_N2C);
qe0 = qc0*qb0 + qc' * qb;
qe = -qc*qb0 + qc0*qb - skew(qc)*qb;

% Define Attitude Control System Specifications -> 4 ortho-skew Reaction Wheels
a1 = [1/sqrt(3)     1/sqrt(3)      1/sqrt(3)]';
a2 = [1/sqrt(3)     -1/sqrt(3)     1/sqrt(3)]';
a3 = [-1/sqrt(3)    1/sqrt(3)      1/sqrt(3)]';
a4 = [-1/sqrt(3)    -1/sqrt(3)     1/sqrt(3)]';
as = [a1, a2, a3, a4];

Is = 1;     % kg m^2
It = 0.5;   % kg m^2

Ja = buildJa(Is, as);       % would be the A matrix for reaction wheels


% Define Mc
Mc = 0;         % no external torques applied on the chaser

% Compute Commanded Torque
Tc = skew(w)*Jc*w - Mc + Jc*wc_dot - Jc*invA*B*wd - sign_qe0_0*Jc*invA*qe;


% Assign State Derivatives
omegas_dot = -Ja' * inv(Ja * Ja') * Tc;
w_dot = inv(Jc) * (Mc - skew(w)*Jc*w - skew(w)*Ja*omegas - Ja*omegas_dot);
qb0_dot = -0.5 * w' * qb;
qb_dot = -0.5 * skew(w) * qb + 0.5 * qb0 * w;
xb_dot = [qb0_dot; qb_dot];

dY(14:17) = xb_dot;
dY(18:20) = w_dot;
dY(21:24) = omegas_dot;


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

