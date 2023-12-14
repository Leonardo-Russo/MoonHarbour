function[dstate, omega_vect, omega_dot_vect, aP_c, aG_c, a34B_c, u, rho_d_vect,...
    rho_d_dot_vect, rho_d_ddot_vect, f_norm] = NaturalControl(t, state, ppEarthMCI, ppSunMCI, muM, muE, muS, ...
    timespan, ppMoonECI, deltaE, psiM, deltaM, ppOmegadotVect, t0, ppXd, kp, DU, TU)

%---------------------------PRE-ALLOCATION---------------------------------
dstate = zeros(13,1); 

%--------------------------- STATE VARIABLES ------------------------------

x1t = state(1);
x2t = state(2);
x3t = state(3);
x4t = state(4);
x5t = state(5);
x6t = state(6);
MEE_target = state(1:6);
x_r = state(7);
y_r = state(8);
z_r = state(9);
vx_r = state(10);
vy_r = state(11);
vz_r = state(12);
x7 = state(13);

rho_c_vect = [x_r; y_r; z_r];
rho_c = norm(rho_c_vect);
rho_c_dot_vect = [vx_r; vy_r; vz_r];


%% ------------------- INTEGRATION OF MEEs OF DSG -------------------------

%------------------------- GRAVITATIONAL MODEL ----------------------------

eta = EtaPar(x2t, x3t, x6t);
G = GravMod(x1t, x2t, x3t, x4t, x5t, x6t, muM, eta);

%----------------------- PERTURBING ACCELERATION --------------------------

a34B = ThirdFourthBody(MEE_target', t, ppEarthMCI, ppSunMCI, muE, muS);
aG = GravPert(MEE_target, ppMoonECI, t, timespan, muM, deltaE, psiM, deltaM);

aP = a34B + aG;

%-------------------------- MEE DERIVATIVES -------------------------------

dMEE = G*aP + [zeros(5,1); sqrt(muM/x1t^3)*eta^2];




%% ------------------------- RELATIVE MOTION ------------------------------
% Target = DSG
% Chaser = general SC, whose initial state is a perturbation of DSG initial state


%------------------------------ TARGET-------------------------------------

% Target COEs
COE_target = Equin2Class(MEE_target);
a = COE_target(1);
e = COE_target(2);
i = COE_target(3);
Omega = COE_target(4);
omega = COE_target(5);
f = COE_target(6);
theta = omega + f;
p = a*(1-e^2);


% Rotation matrix from MCI to target LVLH
R_MCI2targetLVLH = R3(theta)*R1(i)*R3(Omega);

% Target perturbations in target LVLH frame
a34B_t = ThirdFourthBody_Relative(MEE_target, COE_target, t, ppEarthMCI, ppSunMCI, muE, muS);
aG_t = GravPert(MEE_target, ppMoonECI, t, timespan, muM, deltaE, psiM, deltaM);
aP_t = a34B_t + aG_t;

f_r = aP_t(1);
f_theta = aP_t(2);
f_h = aP_t(3);

% Target state in MCI
targetState = Class2Cart(COE_target,muM);
r_target_MCI = targetState(1:3);
v_target_MCI = targetState(4:6);
r_t = norm(r_target_MCI);

r_target_lvlh = R_MCI2targetLVLH*r_target_MCI;

h_vect = cross(r_target_MCI, v_target_MCI);
h = norm(h_vect);

q = ( rho_c^2 + 2*dot(rho_c_vect,r_target_lvlh) ) / r_t^2 ;

% Time derivative of COEs used in rotation matrix
i_dot = r_t*cos(theta)/h*f_h;
Omega_dot = r_t*sin(theta)/(h*sin(i))*f_h;
omega_dot = -r_t*sin(theta)*cos(i)/(h*sin(i))*f_h - cos(f)/e*sqrt(p/muM)*f_r + sqrt(p/muM)*sin(f)*( e*cos(f) + 2 )/( e*(1+e*cos(f)) )*f_theta;
f_dot = sqrt(muM/p^3)*( 1+e*cos(f) )^2 + cos(f)/e*sqrt(p/muM)*f_r - sqrt(p/muM)*sin(f)*( e*cos(f) + 2 )/( e*(1+e*cos(f)) )*f_theta;
theta_dot = omega_dot + f_dot;

% Angular velocity of LVLH frame as seen from MCI
r_hat = [1;0;0];
theta_hat = [0;1;0];
h_hat = [0;0;1];


omega_vect = ( Omega_dot*sin(i)*sin(theta) + i_dot*cos(theta) )*r_hat + ...
             ( Omega_dot*sin(i)*cos(theta) - i_dot*sin(theta) )*theta_hat + ...
             ( Omega_dot*cos(i) + theta_dot )*h_hat;

% Computing time derivative of omega_vect at current integration time
omega_dot_vect = ppsval(ppOmegadotVect, t);



%------------------------------ CHASER ------------------------------------


% Chaser state in MCI frame
f_pert = aP_t;
[r_chaser_MCI,v_chaser_MCI] = lvlh2inertial(rho_c_vect,rho_c_dot_vect,COE_target,f_pert,muM);
r_c = norm(r_chaser_MCI);

% Chaser COEs
COE_chaser = Cart2Class([r_chaser_MCI;v_chaser_MCI], muM);
a_c = COE_chaser(1);
e_c = COE_chaser(2);
i_c = COE_chaser(3);
Omega_c = COE_chaser(4);
omega_c = COE_chaser(5);
f_c = COE_chaser(6);
theta_c = omega_c + f_c;

% Rotation matrix from chaser LVLH to MCI
R_chaserLVLH2MCI = R3(-Omega_c)*R1(-i_c)*R3(-theta_c);

% Chaser MEEs
MEE_chaser = Class2Equin(COE_chaser);


% Chaser perturbations in target LVLH frame
a34B_c = ThirdFourthBody_Relative(MEE_chaser, COE_target, t, ppEarthMCI, ppSunMCI, muE, muS);
aG = GravPert(MEE_chaser, ppMoonECI, t, timespan, muM, deltaE, psiM, deltaM);
aG_c = R_MCI2targetLVLH*R_chaserLVLH2MCI*aG;

aP_c = a34B_c + aG_c;




%--------------------- FEEDBACK LINEARIZATION -----------------------------
f = -2*cross(omega_vect, rho_c_dot_vect) - cross(omega_dot_vect, rho_c_vect) - cross(omega_vect, cross(omega_vect,rho_c_vect)) + ...
           muM/r_t^3* q*(2+q+(1+q)^0.5)/(  (1+q)^1.5 *( (1+q)^0.5+1 )  )*r_target_lvlh ...
          -muM/r_c^3* rho_c_vect + aP_c - aP_t ;
f_norm = norm(f);

Kp = kp*eye(3,3);
%condition of critical dumping
zeta = 1; 
kd = 2*zeta*sqrt(kp);
Kd = kd*eye(3,3);

rho_d_vect = [ppval(ppXd(1),t)
              ppval(ppXd(2),t)
              ppval(ppXd(3),t)];

pp_rho_dot = [fnder(ppXd(1),1)
              fnder(ppXd(2),1)
              fnder(ppXd(3),1)];
rho_d_dot_vect = [ppval(pp_rho_dot(1),t)
                  ppval(pp_rho_dot(2),t)
                  ppval(pp_rho_dot(3),t)];

pp_rho_ddot = [fnder(pp_rho_dot(1),1)
               fnder(pp_rho_dot(2),1)
               fnder(pp_rho_dot(3),1)];
rho_d_ddot_vect = [ppval(pp_rho_ddot(1),t)
                   ppval(pp_rho_ddot(2),t)
                   ppval(pp_rho_ddot(3),t)];


u = -f + rho_d_ddot_vect - Kd*(rho_c_dot_vect-rho_d_dot_vect) - Kp*(rho_c_vect -rho_d_vect);


%--------------------- DIFFERENTIAL EQUATIONS ----------------------------- 
x_r_dot = vx_r;
y_r_dot = vy_r;
z_r_dot = vz_r;
v_r_dot = f + u;

% Mass ratio diff. eq.
c = 30/DU*TU;           % effective exhaust velocity = 30 km/s
x7_dot = - x7*norm(u)/c;

dstate(:,1) = [dMEE; x_r_dot; y_r_dot; z_r_dot; v_r_dot; x7_dot];

end

