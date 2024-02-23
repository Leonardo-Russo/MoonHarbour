function theta_f = finalAngle(rho_0, rho_f, rho_dot_f, tf, t0, theta_0)

    if abs(rho_0-rho_f) == 0
        % In this case the spline is a series of zeroes, so its derivative
        % is also 0 because it's a flat distribution. For this reason, it
        % doesn't exist a best final angle, closest to the desired final
        % velocity. So an arbitrary final angle, namely 0, is chosen.
        theta_f = 0;
    elseif rho_dot_f == 0
        % For the properties of Chebischev distribution, a final angle
        % equal to 0 guarantees final velocity very close to 0, beacuse it 
        % produces a flat ending distribution of points (with derivative ~ 0). 
        % Then 0 is chosen as final angle.
        theta_f = 0;
    else
        if theta_0 > pi/2
            rho_dot_max = - (rho_f - rho_0)/(cos(pi/2) - cos(theta_0))*sin(pi/2)*(pi/2 - theta_0)/(tf - t0);
        else
            % 
            % rho_dot_max = - (rho_f - rho_0)/(cos(theta_0) - cos(theta_0))*sin(pi/2)*(theta_0 - theta_0)/(tf - t0);
        end
        if abs(rho_dot_max) > abs(rho_dot_f) && sign(rho_dot_max) == sign(rho_dot_f)
            % N = 1e6;
            % theta_f = flip(0:pi/N:pi);
            % len = length(theta_f);
            % rho_dot = zeros(len,1);
            % for k = 1 : len
            %     rho_dot(k) = - (rho_f - rho_0)/(cos(theta_f(k)) - cos(theta_0))*sin(theta_f(k))*(theta_f(k) - theta_0)/(tf - t0);
            % end
            % 
            % % One looks for the closest value of rho_dot to the desired rho_dot_f
            % closest = interp1(rho_dot(round(N/2):end), rho_dot(round(N/2):end), rho_dot_f, 'nearest', 'extrap');
            % 
            % % The index corresponding to this value is retrieved
            % index = find(rho_dot == closest);
            % 
            % theta_f = theta_f(index);

            syms  theta_f
            eqn = rho_dot_f == - (rho_f - rho_0)/(cos(theta_f) - cos(theta_0))*sin(theta_f)*(theta_f - theta_0)/(tf - t0);
            theta_f = double(vpasolve(eqn,theta_f));
        else
            theta_f = pi/2;
        end
    end

end