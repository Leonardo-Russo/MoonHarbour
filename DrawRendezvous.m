function DrawRendezvous(Xt_MCI, Xc_MCI, bookmark, opt)
% Description: Create a 3D animation of the rendezvous procedure and save it as an MP4 video.

videoFilename = 'rendezvous.mp4';
fps = 30; % Frames per second for the video

Xt = Xt_MCI(:, 1);
Yt = Xt_MCI(:, 2);
Zt = Xt_MCI(:, 3);

Xc = Xc_MCI(:, 1);
Yc = Xc_MCI(:, 2);
Zc = Xc_MCI(:, 3);

% Define the Initial Point of the Trajectories for Target and State
TrajT = plot3(Xt(1), Yt(1), Zt(1), 'Color', '#d1d1d1', 'Linestyle', '-.', 'LineWidth', 1.5);
hold on
TrajC = plot3(Xc(1), Yc(1), Zc(1), 'Color', '#ff7403', 'Linestyle', '-', 'LineWidth', 1.5);
T = plot3(Xt(1), Yt(1), Zt(1), 'Color', '#3232a8', 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 15);
C = plot3(Xc(1), Yc(1), Zc(1), 'Color', '#d6491e', 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 15);

grid on
axis equal
xlabel('$x$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$y$', 'interpreter', 'latex', 'fontsize', 12)
zlabel('$z$', 'interpreter', 'latex', 'fontsize', 12)
view(3)

N = size(Xt_MCI, 1);

if opt.saveplots
    % Set up the video writer
    writerObj = VideoWriter(videoFilename, 'MPEG-4');
    writerObj.FrameRate = fps;
    open(writerObj);
end

for i = 1 : N
    
    % Set camera distance
    alpha = 1.2;
    camdist = alpha * norm(Xt_MCI(i, 1:3) - Xc_MCI(i, 1:3));
    
    % Update the Plot
    set(TrajT, 'XData', Xt(1:i), 'YData', Yt(1:i), 'ZData', Zt(1:i));
    set(TrajC, 'XData', Xc(1:i), 'YData', Yc(1:i), 'ZData', Zc(1:i));
    set(T, 'XData', Xt(i), 'YData', Yt(i), 'ZData', Zt(i));
    set(C, 'XData', Xc(i), 'YData', Yc(i), 'ZData', Zc(i));
    
    % Adjust axes limits to follow the TrajT
    xlim([Xc(i) - camdist, Xc(i) + camdist]);
    ylim([Yc(i) - camdist, Yc(i) + camdist]);
    zlim([Zc(i) - camdist, Zc(i) + camdist]);

    if opt.saveplots
        % Capture the plot as a frame and write to the video file
        frame = getframe(gcf);
        writeVideo(writerObj, frame);
    else
        pause(1/fps);
    end

end

if opt.saveplots
    close(writerObj);
end

end

