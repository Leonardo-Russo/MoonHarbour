function [ppXd, ViaPoints, t2] = ChebyschevReferenceTrajectory(initialRelativeState, finalRelativeState, t0, tf)
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

    % Initial conditions
    rho_r0 = initialRelativeState(7);
    rho_theta0 = initialRelativeState(8);
    rho_h0 = initialRelativeState(9);
    rho_dot_r0 = 0;
    rho_dot_theta0 = 0;
    rho_dot_h0 = 0;

    % Final conditions
    rho_rf = finalRelativeState(7); 
    rho_thetaf = finalRelativeState(8);
    rho_hf = finalRelativeState(9);
    rho_dot_rf = finalRelativeState(10);
    rho_dot_thetaf = finalRelativeState(11);
    rho_dot_hf = finalRelativeState(12);

    % Generation of via points
    theta0 = pi;
    finalAngle_r = finalAngle(rho_r0,rho_rf, rho_dot_rf, tf, t0,theta0);
    finalAngle_theta = finalAngle(rho_theta0,rho_thetaf, rho_dot_thetaf,tf,t0,theta0);
    finalAngle_h = finalAngle(rho_h0,rho_hf, rho_dot_hf, tf, t0,theta0);

    % N2 = fix((tf-t0)*TU/180);
    N2 = 2;
    t2 = zeros(N2,1);
    for i = 1 : N2
        % time of each via point
        t2(i) = t0 + (i-1)/(N2-1) * (tf-t0);
    end
    rho_r = Chebspace(rho_r0, rho_rf, pi, finalAngle_r, N2)';
    rho_theta = Chebspace(rho_theta0, rho_thetaf, pi, finalAngle_theta, N2)';
    rho_h = Chebspace(rho_h0, rho_hf, pi, finalAngle_h, N2)';



    % Reference trajectory (spline)
  
    pp_rho_r = csape(t2,[rho_dot_r0 rho_r' rho_dot_rf], [0,1]); 
    pp_rho_theta = csape(t2,[rho_dot_theta0 rho_theta' rho_dot_thetaf], [0,1]); 
    pp_rho_h = csape(t2,[rho_dot_h0 rho_h' rho_dot_hf], [0,1]); 

    ppXd = [pp_rho_r; pp_rho_theta; pp_rho_h];
    ViaPoints = [rho_r, rho_theta, rho_h];

end