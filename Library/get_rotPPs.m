function R_PPs = get_rotPPs(times, R)
% Description: this function interpolates a rotation matrix element by
% element using splines. The rotation matrix must be provided as a 
% 3x3xN element.

R_PPs = [];

for j = 1 : 3
    for i = 1 : 3       
        PP_ij = spline(times, R(i, j, :));
        R_PPs = [R_PPs; PP_ij];
    end
end

R_PPs = reshape(R_PPs, 3, 3);

end