function Rdot_PPs = fnder_rots(R_PPs)
% Description: this function derives the interpolating structure of a 
% rotation matrix element by element. The rotation matrix must be provided 
% as a 3x3 PP structure.

Rdot_PPs = [];

for j = 1 : 3
    for i = 1 : 3       
        PP_ij = fnder(R_PPs(i, j), 1);
        Rdot_PPs = [Rdot_PPs; PP_ij];
    end
end

Rdot_PPs = reshape(Rdot_PPs, 3, 3);

end