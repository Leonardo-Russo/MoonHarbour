function P = DrawTrajMCI3D(rMatrixMCI, color, linestyle)
% Description: create a 3D Plot of the propagated orbit.

if nargin < 2
    color = '#ff7403';
end

if nargin < 3
    linestyle = '-';
end

X = rMatrixMCI(:, 1);
Y = rMatrixMCI(:, 2);
Z = rMatrixMCI(:, 3);

[x,y,z]=sphere;

rM = 1738.1;    % km
I = imread('moon.jpg');
surface(rM*x, rM*y, rM*z, flipud(I), 'FaceColor', 'texturemap', 'EdgeColor', 'none', 'CDataMapping', 'direct')

hold on
P = plot3(X,Y,Z,'Color',color, 'Linestyle', linestyle, 'LineWidth', 1.5);
% P = plot3(X,Y,Z,'Color',color, 'Linestyle', 'none', 'LineWidth', 1.5, 'Marker','.', 'MarkerSize',10);

hold off
grid on
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
view([30, 30])


end