function is_safe = check_min_distance(dist, DU, M_ctrl_DA, M_ctrl, min_dist)
% Description: this is a helper function to verify that the emergency
% sphere is not cut by too much. The tolerance is set to 20cm, meaning that
% if the distance is below 9.8m, then it will flag this.

if nargin < 5
    min_dist = 9.8;
end

% fprintf("'dist' is a %.0fx1 vector.\n\n", length(dist));

dist = dist * DU;

is_safe = 1;
for i = 1 : M_ctrl_DA + M_ctrl - 1

    if dist(i) < min_dist * 1e-3
        is_safe = 0;
        % fprintf('crossed minimum distance at index: %.0f\n', i);
    end

end

end

