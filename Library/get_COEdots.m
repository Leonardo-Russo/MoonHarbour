function COEdots = get_COEdots(COE, X_MCI, ap_LVLH)
% Description: this function computes the derivatives of the COE.
% 
% Inputs:
% COE = Classical Orbital Elements
% X_MCI = Position Velocity State in MCI
% ap = perturbing accelerations in LVLH

% Recall Global Variables
global muM Rm DU TU

% Retrieve Data from Input
a = COE(1);
e = COE(2);
incl = COE(3);
Omega = COE(4);
omega = COE(5);
nu = COE(6);
theta_t = COE(5) + COE(6);

rt_MCI = X_MCI(1:3);
vt_MCI = X_MCI(4:6);

% Introduce Local Variables
p = a*(1-e^2);
r = norm(rt_MCI);
h_vect = cross(rt_MCI, vt_MCI);
h = norm(h_vect);

% Retrieve the Perturbing Accelerations
ap_r = ap_LVLH(1);
ap_t = ap_LVLH(2);
ap_h = ap_LVLH(3);

% Compute COEdots - Pontani's Notes
a_dot = 2*a^2/h*(e*sin(nu)*ap_r + p/r*ap_t);
e_dot = sqrt(p/muM)*ap_r*sin(nu) + sqrt(p/muM)*ap_t*(e + e*cos(nu)^2 + 2*cos(nu))/(1 + e*cos(nu));
incl_dot = r*ap_h*cos(theta_t)/h;
Omega_dot = r*ap_h*sin(theta_t)/(h*sin(incl));
omega_dot = -r*ap_h*(sin(theta_t)*cos(incl))/(h*sin(incl)) - ap_r/e*cos(nu)*sqrt(p/muM) + ap_t*sqrt(p/muM)*sin(nu)*(e*cos(nu) + 2)/(e*(1 + e*cos(nu)));
nu_dot = sqrt(muM/p^3)*(1 + e*cos(nu))^2 + ap_r/e*cos(nu)*sqrt(p/muM) - ap_t*sqrt(p/muM)*sin(nu)*(e*cos(nu) + 2)/(e*(1 + e*cos(nu)));

COEdots = [a_dot, e_dot, incl_dot, Omega_dot, omega_dot, nu_dot];

end