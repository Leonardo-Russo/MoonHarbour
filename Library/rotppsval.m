function R = rotppsval(R_PPs, time)
% Description: this function takes an interpolated rotation matrix and
% retrieves its value at a given epoch.

R = eye(3);

for i = 1 : 3
    for j = 1 : 3

        R(i, j) = ppval(R_PPs(i, j), time);

    end
end


end