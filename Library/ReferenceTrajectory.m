function [RHOrefPPs, viapoints, tspan_viapoints, rho_1, t1, l_hat] = ReferenceTrajectory(TC0, TCf, t0, tf, BC, rho_1, t1, l_hat)
% Description: this function computes reference trajectory that should be followed by
% the chaser during the approach to the target (rendezvous). This is
% carried out in target LVLH frame, so we work with relative distances
% (rho = r_chaser - r_target). 
% The reference trajectory is generated considering a Chebyschev distribution
% for via points of each component of relative position rho, bringing rho(t0)
% components to 0 in 12h. This distribution has the peculiar feature that
% can produce a denser set of points at the beginning (where dynamic is
% faster) and a less dense set of point at the end (where dynamic is
% slower). In this way N via points are generated, which are then 
% connected with a spline.
% The reference trajectory is then evaluated at a test timespan and if it
% crosses the emergency radius, then a new procedure will begin to avoid
% it by setting a via point tangent to the sphere.
%
% t0 and tf must be nondimensional.
% TU must be expressed in seconds.

if nargin < 6
    rho_1 = NaN;
end
if nargin < 7
    t1 = NaN;
end
if nargin < 8
    l_hat = NaN;
end

global DU

% Initial conditions
rho_r0 = TC0(7);
rho_t0 = TC0(8);
rho_h0 = TC0(9);
rhodot_r0 = TC0(10);
rhodot_t0 = TC0(11);
rhodot_h0 = TC0(12);
rho_0 = TC0(7:9);
rhodot_0 = TC0(10:12);

% Final conditions
rho_rf = TCf(7); 
rho_tf = TCf(8);
rho_hf = TCf(9);
rhodot_rf = TCf(10);
rhodot_tf = TCf(11);
rhodot_hf = TCf(12);
rho_f = TCf(7:9);
rhodot_f = TCf(10:12);

% Generation of via points
theta0 = pi;
finalAngle_r = finalAngle(rho_r0, rho_rf, rhodot_rf, tf, t0, theta0);
finalAngle_t = finalAngle(rho_t0, rho_tf, rhodot_tf, tf, t0, theta0);
finalAngle_h = finalAngle(rho_h0, rho_hf, rhodot_hf, tf, t0, theta0);

% Free Unconstrained Trajectory
N = 2;          % n° of viapoints

tspan_viapoints = zeros(N,1);
for i = 1 : N
    tspan_viapoints(i) = t0 + (i-1)/(N-1) * (tf-t0);    % epoch of each via point
end

rho_r = Chebspace(rho_r0, rho_rf, pi, finalAngle_r, N)';
rho_t = Chebspace(rho_t0, rho_tf, pi, finalAngle_t, N)';
rho_h = Chebspace(rho_h0, rho_hf, pi, finalAngle_h, N)';

% Compute Reference trajectory (spline)
rho_rPPs = csape(tspan_viapoints,[rhodot_r0 rho_r' rhodot_rf], BC); 
rho_tPPs = csape(tspan_viapoints,[rhodot_t0 rho_t' rhodot_tf], BC); 
rho_hPPs = csape(tspan_viapoints,[rhodot_h0 rho_h' rhodot_hf], BC); 

RHOrefPPs = [rho_rPPs; rho_tPPs; rho_hPPs];
viapoints = [rho_r, rho_t, rho_h];


M = 1000;   % n° of points for sample tspan
tspan_check = linspace(t0, tf, M)';
dist = zeros(M, 1);

sphere_radius = 10e-3/DU;           % 12m of emergency sphere radius
emergency_radius = 10e-3/DU;

% Check for Emergency Sphere Intersection
collision = 0;    
for j = 1 : M
    dist(j) = norm(ppsval(RHOrefPPs, tspan_check(j)));
    if dist(j) <= emergency_radius
        collision = 1;
    end
end


% Three Via Points Method
if collision

    % fprintf('Switching to 3 Via Points Method.\n')

    % Compute MidTime
    if isnan(t1)
        [~, min_idx] = min(dist);
        t1 = tspan_check(min_idx);
    end
    
    % Compute Auxiliary Reference Frame
    if isnan(l_hat)
        l_hat = cross(rho_0, rho_f) / norm(cross(rho_0, rho_f));
    end

    lambda0_hat = cross(l_hat, rho_0)/norm(rho_0);
    rho_0_hat = rho_0 / norm(rho_0);
    xi = acos(sphere_radius/norm(rho_0));

    % Compute rho_1
    if isnan(rho_1)
        rho_1 = sphere_radius * (cos(xi)*rho_0_hat + sin(xi)*lambda0_hat);
    end

    % Create Additional Via Point
    rho_r1 = rho_1(1);
    rho_t1 = rho_1(2);
    rho_h1 = rho_1(3);
    
    rho_r = [rho_r0, rho_r1, rho_rf]';
    rho_t = [rho_t0, rho_t1, rho_tf]';
    rho_h = [rho_h0, rho_h1, rho_hf]';

    if t1 > t0
    
        % Compute Reference Trajectory
        TangentPPs = TangentInterpolation(rho_0, rho_f, rho_1, rhodot_0, rhodot_f, t0, tf, t1, l_hat);
    
        if ~isempty(TangentPPs)

            RHOrefPPs = TangentPPs;
            viapoints = [rho_r, rho_t, rho_h];
            tspan_viapoints = [t0, t1, tf]';

        else    % this is Pre-Emergency Manoeuvre
            
            viapoints = [rho_r, rho_t, rho_h];
            tspan_viapoints = [t0, t1, tf]';

            rho_rPPs = csape(tspan_viapoints,[rhodot_r0 rho_r' rhodot_rf], [0, 1]); 
            rho_tPPs = csape(tspan_viapoints,[rhodot_t0 rho_t' rhodot_tf], [0, 1]); 
            rho_hPPs = csape(tspan_viapoints,[rhodot_h0 rho_h' rhodot_hf], [0, 1]); 
            
            RHOrefPPs = [rho_rPPs; rho_tPPs; rho_hPPs];

        end

    end

end


% Compute the Reference Trajectory Velocity and Acceleration
rhodot_rPPs = fnder(RHOrefPPs(1), 1);
rhodot_tPPs = fnder(RHOrefPPs(2), 1);
rhodot_hPPs = fnder(RHOrefPPs(3), 1);
rhoddot_rPPs = fnder(rhodot_rPPs, 1);
rhoddot_tPPs = fnder(rhodot_tPPs, 1);
rhoddot_hPPs = fnder(rhodot_hPPs, 1);

RHOrefPPs = [RHOrefPPs; rhodot_rPPs; rhodot_tPPs; rhodot_hPPs; rhoddot_rPPs; rhoddot_tPPs; rhoddot_hPPs];


end