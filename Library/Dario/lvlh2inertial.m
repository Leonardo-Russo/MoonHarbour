function [r_chaser_inertial, v_chaser_inertial] = lvlh2inertial(r_chaser_lvlh,v_chaser_lvlh, COE_target, f_pert, mu)
% HELP
%
% This function transforms inertial position adn velocity vectors from an
% inertial RF to a LVLH reference frame, under the assumption that the SC
% on which the LVLH frame is a real SC, whose orbit is affected by 
% perturbations.
%
% INPUTS
%   - r_inertial: inertial position vector to be transformed
%   - v_inertial: inertial velocity vector to be transformed
%   - COE: Classical Orbital Elements of the SC on which the LVLH frame is centered
%   - f_pert: perturbing acceleration, expressed in LVLH frame
%   - mu: gravitational parameter of the central body
%
% OUTPUTS
%   -r_lvlh: position vector in LVLH frame
%   -v_lvlh: velocity vector in LVLH frame

a = COE_target(1);         % semi-major axis [km]
e = COE_target(2);         % eccentricity
i = COE_target(3);         % inclination [rad]
Omega = COE_target(4);     % RAAN [rad]
omega = COE_target(5);     % argument of pericenter [rad]
f = COE_target(6);         % true anomaly [rad]
theta = omega + f;  % argument of latitude [rad]

f_r = f_pert(1);
f_theta = f_pert(2);
f_h = f_pert(3);


% Retrieving position an velocity of the LVLH frame (coincident with
% position and velocity of the SC on which the LVLH frame is centered)
inertialTargetState = Class2Cart(COE_target, mu);
r_target_inertial = inertialTargetState(1:3);
v_target_inertial = inertialTargetState(4:6);

r = norm(r_target_inertial);
h_vect = cross(r_target_inertial, v_target_inertial);
h = norm(h_vect);



% Rotation matrix from LVLH frame to inertial frame
% R = R3(-Omega) * R1(-i) * R3(-omega - f) 
% R = [   cos(omega + f)*cos(Omega) - sin(omega + f)*sin(Omega)*cos(i), - sin(omega + f)*cos(Omega) - cos(omega + f)*sin(Omega)*cos(i),  sin(Omega)*sin(i)
%         cos(omega + f)*sin(Omega) + sin(omega + f)*cos(Omega)*cos(i),   cos(omega + f)*cos(Omega)*cos(i) - sin(omega + f)*sin(Omega), -cos(Omega)*sin(i)
%                                           sin(omega + f)*sin(i),                                             cos(omega + f)*sin(i),             cos(i)];
%
%R = R3(-Omega)*R1(-i)*R3(-theta)
R = [cos(Omega)*cos(theta) - sin(Omega)*cos(i)*sin(theta), - cos(Omega)*sin(theta) - sin(Omega)*cos(i)*cos(theta),  sin(Omega)*sin(i)
    sin(Omega)*cos(theta) + cos(Omega)*cos(i)*sin(theta),   cos(Omega)*cos(i)*cos(theta) - sin(Omega)*sin(theta), -cos(Omega)*sin(i)
                                            sin(i)*sin(theta),                                               cos(theta)*sin(i),                cos(i)];


% Time derivative of true anomaly
p = a*(1-e^2);
i_dot = r*cos(theta)/h*f_h;
Omega_dot = r*sin(theta)/(h*sin(i))*f_h;
omega_dot = -r*sin(theta)*cos(i)/(h*sin(i))*f_h - cos(f)/e*sqrt(p/mu)*f_r + sqrt(p/mu)*sin(f)*( e*cos(f) + 2 )/( e*(1+e*cos(f)) )*f_theta;
f_dot = sqrt(mu/p^3)*( 1+e*cos(f) )^2 + cos(f)/e*sqrt(p/mu)*f_r - sqrt(p/mu)*sin(f)*( e*cos(f) + 2 )/( e*(1+e*cos(f)) )*f_theta;
theta_dot = omega_dot + f_dot;


% Time derivative of the rotation matrix
% R_dot = [sin(Omega)*sin(i)*sin(f + omega)*i_dot - cos(Omega)*sin(f + omega)*(f_dot + omega_dot) - sin(Omega)*cos(i)*cos(f + omega)*(f_dot + omega_dot) - sin(Omega)*cos(f + omega)*Omega_dot - cos(Omega)*cos(i)*sin(f + omega)*Omega_dot,            sin(Omega)*sin(f + omega)*Omega_dot - cos(Omega)*cos(f + omega)*(f_dot + omega_dot) + sin(Omega)*sin(i)*cos(f + omega)*i_dot + sin(Omega)*cos(i)*sin(f + omega)*(f_dot + omega_dot) - cos(Omega)*cos(i)*cos(f + omega)*Omega_dot,             cos(Omega)*sin(i)*Omega_dot + sin(Omega)*cos(i)*i_dot
%         cos(Omega)*cos(f + omega)*Omega_dot - sin(Omega)*sin(f + omega)*(f_dot + omega_dot) - sin(Omega)*cos(i)*sin(f + omega)*Omega_dot + cos(Omega)*cos(i)*cos(f + omega)*(f_dot + omega_dot) - cos(Omega)*sin(i)*sin(f + omega)*i_dot ,          - cos(Omega)*sin(f + omega)*Omega_dot - sin(Omega)*cos(f + omega)*(f_dot + omega_dot) - cos(Omega)*sin(i)*cos(f + omega)*i_dot - cos(Omega)*cos(i)*sin(f + omega)*(f_dot + omega_dot) - sin(Omega)*cos(i)*cos(f + omega)*Omega_dot,             sin(Omega)*sin(i)*Omega_dot - cos(Omega)*cos(i)*i_dot
%                                                                                                                                                                  cos(i)*sin(f + omega)*i_dot + sin(i)*cos(f + omega)*(f_dot + omega_dot) ,                                                                                                                                                                     cos(i)*cos(f + omega)*i_dot - sin(i)*sin(f + omega)*(f_dot + omega_dot),                                                     -sin(i)*i_dot];

R_dot =[sin(Omega)*sin(i)*sin(theta)*i_dot - cos(Omega)*sin(theta)*theta_dot - cos(Omega)*cos(i)*sin(theta)*Omega_dot - sin(Omega)*cos(i)*cos(theta)*theta_dot - sin(Omega)*cos(theta)*Omega_dot,   sin(Omega)*sin(theta)*Omega_dot - cos(Omega)*cos(theta)*theta_dot - cos(Omega)*cos(i)*cos(theta)*Omega_dot + sin(Omega)*cos(theta)*sin(i)*i_dot + sin(Omega)*cos(i)*sin(theta)*theta_dot, cos(Omega)*sin(i)*Omega_dot + sin(Omega)*cos(i)*i_dot
        cos(Omega)*cos(theta)*Omega_dot - sin(Omega)*sin(theta)*theta_dot - sin(Omega)*cos(i)*sin(theta)*Omega_dot + cos(Omega)*cos(i)*cos(theta)*theta_dot - cos(Omega)*sin(i)*sin(theta)*i_dot, - cos(Omega)*sin(theta)*Omega_dot - sin(Omega)*cos(theta)*theta_dot - sin(Omega)*cos(i)*cos(theta)*Omega_dot - cos(Omega)*cos(theta)*sin(i)*i_dot - cos(Omega)*cos(i)*sin(theta)*theta_dot, sin(Omega)*sin(i)*Omega_dot - cos(Omega)*cos(i)*i_dot
                                                                                                                                                                                      cos(i)*sin(theta)*i_dot + cos(theta)*sin(i)*theta_dot,                                                                                                                                                                                         cos(i)*cos(theta)*i_dot - sin(i)*sin(theta)*theta_dot,                                                          -sin(i)*i_dot];
 

r_chaser_inertial = R*r_chaser_lvlh + r_target_inertial;
v_chaser_inertial = R*v_chaser_lvlh + R_dot*r_chaser_lvlh + v_target_inertial;



end