function RHOrefPPs = TangentInterpolation(rho_0, rho_f, rho_1, rhodot_0, rhodot_f, t0, tf, t1, l_hat)
% Description: this function performs the custom interpolation following
% our requirements.

% Define the symbolic variables
syms a0_1 a1_1 a2_1 a3_1 a0_2 a1_2 a2_2 a3_2 t real
syms b0_1 b1_1 b2_1 b3_1 b0_2 b1_2 b2_2 b3_2 real
syms c0_1 c1_1 c2_1 c3_1 c0_2 c1_2 c2_2 c3_2 real
vars = [a0_1, a1_1, a2_1, a3_1, a0_2, a1_2, a2_2, a3_2, b0_1, b1_1, b2_1, b3_1, b0_2, b1_2, b2_2, b3_2, c0_1, c1_1, c2_1, c3_1, c0_2, c1_2, c2_2, c3_2];

% Solution Existance Flag
exists_solution = 1;

% Define the Spline Equations
x0_eq = a0_1 + a1_1 * (t - t0) + a2_1 * (t - t0)^2 + a3_1 * (t - t0)^3;
u0_eq = a1_1 + 2 * a2_1 * (t - t0) + 3 * a3_1 * (t - t0)^2;
x1_eq = a0_2 + a1_2 * (t - t1) + a2_2 * (t - t1)^2 + a3_2 * (t - t1)^3;
u1_eq = a1_2 + 2 * a2_2 * (t - t1) + 3 * a3_2 * (t - t1)^2;
ax0_eq = 2*a2_1 + 6*a3_1*(t - t0);
ax1_eq = 2*a2_2 + 6*a3_2*(t - t1);

y0_eq = b0_1 + b1_1 * (t - t0) + b2_1 * (t - t0)^2 + b3_1 * (t - t0)^3;
v0_eq = b1_1 + 2 * b2_1 * (t - t0) + 3 * b3_1 * (t - t0)^2;
y1_eq = b0_2 + b1_2 * (t - t1) + b2_2 * (t - t1)^2 + b3_2 * (t - t1)^3;
v1_eq = b1_2 + 2 * b2_2 * (t - t1) + 3 * b3_2 * (t - t1)^2;
ay0_eq = 2*b2_1 + 6*b3_1*(t - t0);
ay1_eq = 2*b2_2 + 6*b3_2*(t - t1);

z0_eq = c0_1 + c1_1 * (t - t0) + c2_1 * (t - t0)^2 + c3_1 * (t - t0)^3;
w0_eq = c1_1 + 2 * c2_1 * (t - t0) + 3 * c3_1 * (t - t0)^2;
z1_eq = c0_2 + c1_2 * (t - t1) + c2_2 * (t - t1)^2 + c3_2 * (t - t1)^3;
w1_eq = c1_2 + 2 * c2_2 * (t - t1) + 3 * c3_2 * (t - t1)^2;
az0_eq = 2*c2_1 + 6*c3_1*(t - t0);
az1_eq = 2*c2_2 + 6*c3_2*(t - t1);

% Recall the Known Variables
x0 = rho_0(1);
y0 = rho_0(2);
z0 = rho_0(3);
u0 = rhodot_0(1);
v0 = rhodot_0(2);
w0 = rhodot_0(3);
x1 = rho_1(1);
y1 = rho_1(2);
z1 = rho_1(3);
xf = rho_f(1);
yf = rho_f(2);
zf = rho_f(3);
uf = rhodot_f(1);
vf = rhodot_f(2);
wf = rhodot_f(3);

% Compute Local Variables
v1_hat = [subs(u1_eq, t, t1), subs(v1_eq, t, t1), subs(w1_eq, t, t1)]/norm([subs(u1_eq, t, t1), subs(v1_eq, t, t1), subs(w1_eq, t, t1)]);
a0 = [subs(ax0_eq, t, t1), subs(ay0_eq, t, t1), subs(az0_eq, t, t1)];
a1 = [subs(ax1_eq, t, t1), subs(ay1_eq, t, t1), subs(az1_eq, t, t1)];
a0_v1 = dot(a0, v1_hat);
a1_v1 = dot(a1, v1_hat);

% Define the Boundary Conditions
equations = [subs(x0_eq, t, t0) == x0, ...     % extremal conditions
             subs(u0_eq, t, t0) == u0, ...
             subs(x1_eq, t, tf) == xf, ...
             subs(u1_eq, t, tf) == uf, ...
             subs(y0_eq, t, t0) == y0, ...
             subs(v0_eq, t, t0) == v0, ...
             subs(y1_eq, t, tf) == yf, ...
             subs(v1_eq, t, tf) == vf, ...
             subs(z0_eq, t, t0) == z0, ...
             subs(w0_eq, t, t0) == w0, ...
             subs(z1_eq, t, tf) == zf, ...
             subs(w1_eq, t, tf) == wf, ...
             subs(x0_eq, t, t1) == x1, ...     % continuity at rho1
             subs(x1_eq, t, t1) == x1, ...
             subs(y0_eq, t, t1) == y1, ...
             subs(y1_eq, t, t1) == y1, ...
             subs(z0_eq, t, t1) == z1, ...
             subs(z1_eq, t, t1) == z1, ...
             subs(u0_eq, t, t1) == subs(u1_eq, t, t1), ...     % continuity at rhodot1
             subs(v0_eq, t, t1) == subs(v1_eq, t, t1), ...
             subs(w0_eq, t, t1) == subs(w1_eq, t, t1), ...
             dot([subs(u1_eq, t, t1), subs(v1_eq, t, t1), subs(w1_eq, t, t1)], rho_1) == 0, ...    % final conditions
             dot([subs(u1_eq, t, t1), subs(v1_eq, t, t1), subs(w1_eq, t, t1)], l_hat) == 0, ...
             a0_v1 == a1_v1];

% Compute the Solution
solution = solve(equations, vars);

% Extract Coefficients for Each Segment
% Assuming solution contains the coefficients in the order you defined in 'vars'
coeffs_x_1 = [double(solution.a3_1), double(solution.a2_1), double(solution.a1_1), double(solution.a0_1)];
coeffs_x_2 = [double(solution.a3_2), double(solution.a2_2), double(solution.a1_2), double(solution.a0_2)];
coeffs_y_1 = [double(solution.b3_1), double(solution.b2_1), double(solution.b1_1), double(solution.b0_1)];
coeffs_y_2 = [double(solution.b3_2), double(solution.b2_2), double(solution.b1_2), double(solution.b0_2)];
coeffs_z_1 = [double(solution.c3_1), double(solution.c2_1), double(solution.c1_1), double(solution.c0_1)];
coeffs_z_2 = [double(solution.c3_2), double(solution.c2_2), double(solution.c1_2), double(solution.c0_2)];

% Check the Existance of a Solution
if isempty(coeffs_x_1) || isempty(coeffs_x_2) || isempty(coeffs_y_1) || isempty(coeffs_y_2) || isempty(coeffs_z_1) || isempty(coeffs_z_2)
    exists_solution = 0;
    warning('No possible solution for the Tangent Interpolation! Proceeding with the two-point interpolation instead.')
end

if exists_solution

    % Breakpoints
    breaks = [t0, t1, tf];
    
    % Create Spline Structure Using mkpp
    rho_rPPs = mkpp(breaks, [coeffs_x_1; coeffs_x_2]);
    rho_tPPs = mkpp(breaks, [coeffs_y_1; coeffs_y_2]);
    rho_hPPs = mkpp(breaks, [coeffs_z_1; coeffs_z_2]);
    RHOrefPPs = [rho_rPPs; rho_tPPs; rho_hPPs];

else

    RHOrefPPs = [];

end


end