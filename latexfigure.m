%% Rotations - Leonardo Russo

close all
clear


addpath('Library/')

%% Core

% Define the eigenaxis and eigenangle
eigenaxis = [1, 1, 1];  % Example axis, should be a unit vector
eigenangle = deg2rad(45);  % Example angle in radians

% Normalize the eigenaxis
eigenaxis = eigenaxis / norm(eigenaxis);

% Quaternion components from eigenaxis and eigenangle
q0 = cos(eigenangle / 2);
q1 = eigenaxis(1) * sin(eigenangle / 2);
q2 = eigenaxis(2) * sin(eigenangle / 2);
q3 = eigenaxis(3) * sin(eigenangle / 2);

R = q2C(q0, [q1, q2, q3]');
R = R1(pi/2) * R;
[q0, q] = C2q(R);
q1 = q(1);
q2 = q(2);
q3 = q(3);

% Quaternion to 313 Euler angles (ZXZ convention)
% Calculate the 313 Euler angles
phi = atan2(2*(q1*q3 + q0*q2), q0^2 + q1^2 - q2^2 - q3^2);  % First rotation around Z
theta = acos(q0^2 - q1^2 - q2^2 + q3^2);  % Rotation around X
psi = atan2(2*(q2*q3 + q0*q1), q0^2 - q1^2 + q2^2 - q3^2);  % Second rotation around Z

% Convert angles from radians to degrees
phi_deg = rad2deg(phi);
theta_deg = rad2deg(theta);
psi_deg = rad2deg(psi);

% Display the results
fprintf('313 Euler angles (degrees):\n');
fprintf('phi (first rotation around Z): %.2f\n', phi_deg);
fprintf('theta (rotation around X): %.2f\n', theta_deg);
fprintf('psi (second rotation around Z): %.2f\n', psi_deg);

% Quaternion to Yaw-Pitch-Roll (ZYX convention)
yaw = atan2(2*(q0*q3 + q1*q2), 1 - 2*(q2^2 + q3^2));  % Rotation around Z
pitch = asin(2*(q0*q2 - q3*q1));  % Rotation around Y
roll = atan2(2*(q0*q1 + q2*q3), 1 - 2*(q1^2 + q2^2));  % Rotation around X

% Convert angles from radians to degrees
yaw_deg = rad2deg(yaw);
pitch_deg = rad2deg(pitch);
roll_deg = rad2deg(roll);

% Display the results
fprintf('\nYaw-Pitch-Roll (degrees):\n');
fprintf('Yaw (rotation around Z): %.2f\n', yaw_deg);
fprintf('Pitch (rotation around Y): %.2f\n', pitch_deg);
fprintf('Roll (rotation around X): %.2f\n', roll_deg);

