function [RHOrefPPs, viapoints, tspan_viapoints] = ChebyschevReferenceTrajectory(TC0, TCf, t0, tf)
% Description: this function computes reference trajectory that should be followed by
% the chaser during the approach to the target (rendezvous). This is
% carried out in target LVLH frame, so we work with relative distances
% (rho = r_chaser - r_target). 
% The reference trajectory is generated considering a Chebyschev distribution
% for via points of each component of relative position rho, bringing rho(t0)
% components to 0 in 12h. This distribution has the peculiar feature that
% can produce a denser set of points at the beginning (where dynamic is
% faster) and a less dense set of point at the end (where dynamic is
% slower).
% In this way N via points are generated, which are then 
% connected with a spline.
%
% t0 and tf must be nondimensional.
% TU must be expressed in seconds.
% 
% Credits: This function was developed building upon a core reference 
% provided by Dario Sanna.
% 
% Note: I have added the trajectory verification step.

global DU

% Initial conditions
rho_r0 = TC0(7);
rho_theta0 = TC0(8);
rho_h0 = TC0(9);
rho_dot_r0 = 0;
rho_dot_theta0 = 0;
rho_dot_h0 = 0;

% Final conditions
rho_rf = TCf(7); 
rho_thetaf = TCf(8);
rho_hf = TCf(9);
rho_dot_rf = TCf(10);
rho_dot_thetaf = TCf(11);
rho_dot_hf = TCf(12);

% Generation of via points
theta0 = pi;
finalAngle_r = finalAngle(rho_r0,rho_rf, rho_dot_rf, tf, t0,theta0);
finalAngle_theta = finalAngle(rho_theta0,rho_thetaf, rho_dot_thetaf,tf,t0,theta0);
finalAngle_h = finalAngle(rho_h0,rho_hf, rho_dot_hf, tf, t0,theta0);

% Check Free Trajectory
N = 2;          % n° of viapoints

tspan_viapoints = zeros(N,1);
for i = 1 : N
    tspan_viapoints(i) = t0 + (i-1)/(N-1) * (tf-t0);    % epoch of each via point
end

rho_r = Chebspace(rho_r0, rho_rf, pi, finalAngle_r, N)';
rho_theta = Chebspace(rho_theta0, rho_thetaf, pi, finalAngle_theta, N)';
rho_h = Chebspace(rho_h0, rho_hf, pi, finalAngle_h, N)';

% Compute Reference trajectory (spline)
pp_rho_r = csape(tspan_viapoints,[rho_dot_r0 rho_r' rho_dot_rf], [0,1]); 
pp_rho_theta = csape(tspan_viapoints,[rho_dot_theta0 rho_theta' rho_dot_thetaf], [0,1]); 
pp_rho_h = csape(tspan_viapoints,[rho_dot_h0 rho_h' rho_dot_hf], [0,1]); 

RHOrefPPs = [pp_rho_r; pp_rho_theta; pp_rho_h];
viapoints = [rho_r, rho_theta, rho_h];


M = 1000;   % n° of points for sample tspan
tspan_check = linspace(t0, tf, M)';

tol = 10e-3/DU;      % 10m of emergency sphere tolerance
fail = 0;

for j = 1 : M
    RHOref_check = ppsval(RHOrefPPs, tspan_check(j));
    if norm(RHOref_check(1:3)) <= tol
        fprintf('Switching to 3 Via Points Method.\n')
        fail = 1;
        break
    end
end

if fail

    % Set Additional Via Point
    N = 3;          % n° of viapoints
    
    tspan_viapoints = zeros(N,1);
    for i = 1 : N
        tspan_viapoints(i) = t0 + (i-1)/(N-1) * (tf-t0);    % epoch of each via point
    end
    
    rho_r = Chebspace(rho_r0, rho_rf, pi, finalAngle_r, N)';
    rho_theta = Chebspace(rho_theta0, rho_thetaf, pi, finalAngle_theta, N)';
    rho_h = Chebspace(rho_h0, rho_hf, pi, finalAngle_h, N)';
    
    % Compute Reference trajectory (spline)
    pp_rho_r = csape(tspan_viapoints,[rho_dot_r0 rho_r' rho_dot_rf], [0,1]); 
    pp_rho_theta = csape(tspan_viapoints,[rho_dot_theta0 rho_theta' rho_dot_thetaf], [0,1]); 
    pp_rho_h = csape(tspan_viapoints,[rho_dot_h0 rho_h' rho_dot_hf], [0,1]); 
    
    RHOrefPPs = [pp_rho_r; pp_rho_theta; pp_rho_h];
    viapoints = [rho_r, rho_theta, rho_h];
    
    
    M = 1000;   % n° of points for sample tspan
    tspan_check = linspace(t0, tf, M)';
    
    tol = 10e-3/DU;      % 10m of emergency sphere tolerance
    
    for j = 1 : M

        RHOref_check = ppsval(RHOrefPPs, tspan_check(j));
        
        if norm(RHOref_check(1:3)) >= tol
            fprintf('Damn.\n')
            break
        end

    end

end


end