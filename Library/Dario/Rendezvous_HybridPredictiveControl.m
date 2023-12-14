function [dstate, omega_vect, omega_dot_vect, aP_c, aG_c, a34B_c, u, rho_d_vect,...
    rho_d_dot_vect, rho_d_ddot_vect, f, f_norm, kp_out, k_type] = Rendezvous_HybridPredictiveControl(t, state, ppEarthMCI, ppSunMCI, muM, muE, muS, ...
    timespan, ppMoonECI, deltaE, psiM, deltaM, ppOmegaVect, t0, ppXd, DU, TU, checkTimes, Delta_t, N_inner_integration, P)

    persistent kp verifiedTimes index stopPrediction u_max

    if t == t0
        verifiedTimes = zeros(length(checkTimes),1);
        index = 1;
        stopPrediction = 0;
        u_max = 0;
    end
    
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
    rho_r = state(7);
    rho_theta = state(8);
    rho_h = state(9);
    rho_dot_r = state(10);
    rho_dot_theta = state(11);
    rho_dot_h = state(12);
    x7 = state(13);
    
    rho_c_vect = [rho_r; rho_theta; rho_h];
    rho_c = norm(rho_c_vect);
    rho_c_dot_vect = [rho_dot_r; rho_dot_theta; rho_dot_h];
    rho_c_dot = norm(rho_c_dot_vect);
    
    
    %% ------------------- INTEGRATION OF MEEs OF DSG -------------------------
    
    %------------------------- GRAVITATIONAL MODEL ----------------------------
    
    eta = EtaPar(x2t, x3t, x6t);
    G = GravMod(x1t, x2t, x3t, x4t, x5t, x6t, muM, eta);
    
    %----------------------- PERTURBING ACCELERATION --------------------------
    
    a34B = ThirdBody(MEE_target, t, ppEarthMCI, ppSunMCI, muE, muS, timespan, muM);
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
    a34B_t = ThirdBody_rel(MEE_target, COE_target, t, ppEarthMCI, ppSunMCI, muE, muS, timespan, muM);
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
    pp_der = [fnder(ppOmegaVect(1), 1)
              fnder(ppOmegaVect(2), 1)
              fnder(ppOmegaVect(3), 1)];
    omega_dot_vect = [ppval(pp_der(1),t)
                      ppval(pp_der(2),t)
                      ppval(pp_der(3),t)];
    
    
    
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
    a34B_c = ThirdBody_rel(MEE_chaser, COE_target, t, ppEarthMCI, ppSunMCI, muE, muS, timespan, muM);
    aG = GravPert(MEE_chaser, ppMoonECI, t, timespan, muM, deltaE, psiM, deltaM);
    aG_c = R_MCI2targetLVLH*R_chaserLVLH2MCI*aG;
    
    aP_c = a34B_c + aG_c;
    
    
    
    
    %--------------------- FEEDBACK LINEARIZATION -----------------------------
    f = -2*cross(omega_vect, rho_c_dot_vect) - cross(omega_dot_vect, rho_c_vect) - cross(omega_vect, cross(omega_vect,rho_c_vect)) + ...
               muM/r_t^3* q*(2+q+(1+q)^0.5)/(  (1+q)^1.5 *( (1+q)^0.5+1 )  )*r_target_lvlh ...
              -muM/r_c^3* rho_c_vect + aP_c - aP_t ;
    f_norm = norm(f);
    
    
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
    
    
    %condition of critical dumping
    zeta = 1; 
    
    g0 = 9.80665;
    u_limit =  5e-5*g0/(1000*DU)*TU^2;
    threshold = 0.9*u_limit;
    
    % Predictive propagation
    if t > checkTimes(index) && verifiedTimes(index) == 0 && stopPrediction == 0

        % Forward integration
       
        t0_int = t;
        tf_int = t + Delta_t;
        [tspan, ControlledRelativeState] = ode_Ham(@(t, state) Rendezvous(t, state, ppEarthMCI, ppSunMCI, muM, ...
             muE, muS, timespan, ppMoonECI, deltaE, psiM, deltaM, ppOmegaVect, t0, ppXd, kp, DU, TU), ...
             [t0_int, tf_int], state, N_inner_integration);

        len = length(tspan);
        u_future = zeros(3,len);
        u_norms = zeros(len,1);

        % Post processing
        for s = 1 : len
            [~, ~, ~, ~, ~, ~,u_future(:,s), ...
            ~, ~, ~, ~] = Rendezvous(tspan(s), ...
            ControlledRelativeState(s,:), ppEarthMCI, ppSunMCI, muM, ...
            muE, muS, timespan, ppMoonECI, deltaE, psiM, deltaM, ppOmegaVect, t0, ppXd, kp, DU, TU);
            u_norms(s) = norm(u_future(:,s));
        end
        
        % if (t-t0)*TU/3600 > 1.2
        %     % Debugging plot
            % figure()
            % Day = 86400;
            % u_limits = u_limit*ones(len,1);
            % p1 = plot((tspan - t0) * TU / Day*24, u_norms*1000*DU/TU^2, 'Color', 'b', 'LineWidth', 1.2);
            % hold on
            % p2 = plot((tspan - t0) * TU / Day*24, u_limits*1000*DU/TU^2, 'r--', 'LineWidth', 1.2);
            % xlabel('t [Hours]')
            % ylabel('$u \left[\frac{m}{s^2}\right]$', 'Interpreter', 'latex')
            % set(gca, 'FontSize', 20)
            % legend([p1,p2], '$Actual \ u$', '$u_{max}$','Location', 'best', 'Fontsize', 12, 'Interpreter','latex');
        % end

        [crosses, indexes, ~] = detectCross(u_norms, threshold);
        if crosses <= 1
            % this is a cross from above, that's ok
            verifiedTimes(index) = 1; 
            index = index +1;
        elseif crosses >= 2
            % this means we have at least a cross from below, i. e. the
            % threshold has been overcome by the feedback reaction 
            stopPrediction = 1;
            stopSaturationTime = t;
            save(P, "stopSaturationTime")
        end
      
        if isempty(indexes)
            ind = 1;
        else
            ind = indexes(1);
        end

        [u_max, ind2] = max(u_norms(ind+1:end)); % first values are excluded, because it's in saturation ( therefore u_norms(1) == u_limit )
        
    end

    if u_max < threshold  
        a = -f + rho_d_ddot_vect;
        b = 2*(rho_c_dot_vect - rho_d_dot_vect);
        c = (rho_c_vect - rho_d_vect);

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

    % end

    else
        % Threshold overcome
        k_type = NaN;
    
    end
    
    kp_out = kp;

    Kp = kp*eye(3,3);
    kd = 2*zeta*sqrt(kp);
    Kd = kd*eye(3,3);
    
    u = -f + rho_d_ddot_vect - Kd*(rho_c_dot_vect-rho_d_dot_vect) - Kp*(rho_c_vect -rho_d_vect);
    
    %--------------------- DIFFERENTIAL EQUATIONS -------------------------
    x_r_dot = rho_dot_r;
    y_r_dot = rho_dot_theta;
    z_r_dot = rho_dot_h;
    v_r_dot = f + u;

    % Mass ratio diff. eq.
    c = 30/DU*TU;           % effective exhaust velocity = 30 km/s
    x7_dot = - x7*norm(u)/c;
    
    dstate(:,1) = [dMEE; x_r_dot; y_r_dot; z_r_dot; v_r_dot; x7_dot];
    
    fprintf('Maneuver time [h]: %.4g\n', (t-t0)*TU/3600)
    % [(t-t0)*TU/3600, norm(u)*1000*DU/TU^2, rho_c_vect(1)*1000*DU, rho_c_dot_vect(1)*1000*DU/TU]


end

