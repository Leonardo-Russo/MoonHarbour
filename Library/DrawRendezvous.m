function DrawRendezvous(Xt_MCI, Xc_MCI, RHO_LVLH, bookmark)
% Description: Create a 3D animation of the rendezvous procedure and save it as an MP4 video.

videoFilename = 'Output/rendezvous.mp4';
fps = 60; % Frames per second for the video

az0 = 140;
azf = az0 + 360*1.5;
el = 30;

Xt = Xt_MCI(:, 1);
Yt = Xt_MCI(:, 2);
Zt = Xt_MCI(:, 3);

Xc = Xc_MCI(:, 1);
Yc = Xc_MCI(:, 2);
Zc = Xc_MCI(:, 3);

rho = RHO_LVLH(:, 1);
theta = RHO_LVLH(:, 2);
h = RHO_LVLH(:, 3);

% Sphere radius for the target
sphereRadius = 0.005;       % assume 5m of radius for Lunar Gateway

% Create a sphere for the target
[xs, ys, zs] = sphere;
xs = sphereRadius * xs;
ys = sphereRadius * ys;
zs = sphereRadius * zs;

% Define the Initial Point of the Trajectory in LVLH
subplot(1, 2, 1)
TrajLVLH = plot3(rho(1), theta(1), h(1), 'Color', '#ff7403', 'Linestyle', '-', 'LineWidth', 1.5);
hold on
T_LVLH = plot3(0, 0, 0, 'Color', '#3232a8', 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 15);
C_LVLH = plot3(rho(1), theta(1), h(1), 'Color', '#d6491e', 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 15);

grid on
axis equal
xlabel('$\rho \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$\theta \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
zlabel('$h \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
view([30, el])
title('Chaser in LVLH')

% Define the Initial Point of the Trajectories for Target and State
subplot(1, 2, 2)
TrajT = plot3(Xt(1), Yt(1), Zt(1), 'Color', '#d1d1d1', 'Linestyle', '-.', 'LineWidth', 1.5);
hold on
TrajC = plot3(Xc(1), Yc(1), Zc(1), 'Color', '#ff7403', 'Linestyle', '-', 'LineWidth', 1.5);
T = surf(Xt(1) + xs, Yt(1) + ys, Zt(1) + zs, 'EdgeColor', 'none', 'FaceColor', '#3232a8'); % Sphere for target
C = plot3(Xc(1), Yc(1), Zc(1), 'Color', '#d6491e', 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 15);
title('Chaser and Target in MCI')

axis equal
grid on
xlabel('$x \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$y \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
zlabel('$z \ [km]$', 'interpreter', 'latex', 'fontsize', 12)
view([az0, el])

N = size(Xt_MCI, 1);
az = linspace(az0, azf, N);


% Set up the video writer
writerObj = VideoWriter(videoFilename, 'MPEG-4');
writerObj.FrameRate = fps;
open(writerObj);

for i = 1 : N
    
    % Set camera distance
    alpha = 1.2;
    camdist = alpha * norm(Xt_MCI(i, 1:3) - Xc_MCI(i, 1:3));
    
    % Update the MCI Plot
    set(TrajT, 'XData', Xt(1:i), 'YData', Yt(1:i), 'ZData', Zt(1:i));
    set(TrajC, 'XData', Xc(1:i), 'YData', Yc(1:i), 'ZData', Zc(1:i));
    set(T, 'XData', Xt(i) + xs, 'YData', Yt(i) + ys, 'ZData', Zt(i) + zs);

    set(C, 'XData', Xc(i), 'YData', Yc(i), 'ZData', Zc(i));

    % Update the LVLH Plot
    set(TrajLVLH, 'XData', rho(1:i), 'YData', theta(1:i), 'ZData', h(1:i));
    set(C_LVLH, 'XData', rho(i), 'YData', theta(i), 'ZData', h(i));
    
    % Adjust axes limits to follow the TrajT
    xlims = [Xc(i) - camdist, Xc(i) + camdist];
    ylims = [Yc(i) - camdist, Yc(i) + camdist];
    zlims = [Zc(i) - camdist, Zc(i) + camdist];
    xlim(xlims);
    ylim(ylims);
    zlim(zlims);

    % Set the XTick, YTick, and ZTick properties
    tickIntervals = 5;
    set(gca, 'XTick', linspace(xlims(1), xlims(2), tickIntervals), ...
             'YTick', linspace(ylims(1), ylims(2), tickIntervals), ...
             'ZTick', linspace(zlims(1), zlims(2), tickIntervals));

    view([az(i), el])

    
    % Capture the plot as a frame and write to the video file
    frame = getframe(gcf);
    writeVideo(writerObj, frame);

end


close(writerObj);


end
