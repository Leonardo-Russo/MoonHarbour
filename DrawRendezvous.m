function DrawRendezvous(Xt_MCI, Xc_MCI, bookmark, opt)
% Description: create a 3D animation of the rendezvous procedure.


gifFilename = 'rendezvous.gif';
fps = 1000;

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

        % Capture the plot as a frame
        frame = getframe(gcf);
        [imind, cm] = rgb2ind(frame.cdata, 256);
        
        % Write to the GIF file
        if i == 1
            imwrite(imind, cm, gifFilename, 'gif', 'Loopcount', inf, 'DelayTime', 1/fps);
        else
            imwrite(imind, cm, gifFilename, 'gif', 'WriteMode', 'append', 'DelayTime', 1/fps);
        end

    else

        pause(1/fps);

    end

end

end
