function P = DrawTrajLVLH3D(rMatrixLVLH, color, linestyle, show_spheres)
% Description: create a 3D Plot of the propagated orbit.

if nargin < 2
    color = '#ff7403';
end

if nargin < 3
    linestyle = '-';
end

if nargin < 4
    show_spheres = true;
end

R = rMatrixLVLH(:, 1);
T = rMatrixLVLH(:, 2);
H = rMatrixLVLH(:, 3);

if show_spheres
    [x,y,z]=sphere;
    rT = 5e-3;      % km - S/C approximated as a sphere of 5 meter radius
    % I = imread('titanium.jpg');
    I = imread('black.jpg');
    surface(rT*x, rT*y, rT*z, flipud(I), 'FaceColor', 'texturemap', 'EdgeColor', 'none', 'CDataMapping', 'direct')
    hold on
    radius = 10e-3;
    surface(radius*x, radius*y, radius*z, 'FaceColor', '#ffcf82', 'EdgeColor', 'none', 'CDataMapping', 'direct', 'FaceAlpha', 0.2)
end

P = plot3(R,T,H,'Color',color, 'Linestyle', linestyle, 'LineWidth', 1.5);
hold on
% P = plot3(X,Y,Z,'Color',color, 'Linestyle', 'none', 'LineWidth', 1.5, 'Marker','.', 'MarkerSize',10);

grid on
axis equal
xlabel('$r \ [km]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\theta \ [km]$', 'Interpreter','latex', 'FontSize', 12)
zlabel('$h \ [km]$', 'Interpreter','latex', 'FontSize', 12)
view([30, 30])


end