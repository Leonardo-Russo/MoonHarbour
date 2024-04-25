function A = buildA(Is, a)
% Description: this function creates the A matrix used in the evaluation of
% omegas_dot from the known value of commanded Torque.

A = [];

Nw = size(a, 2);

for i = 1 : Nw
    A = [A, Is*a(:, i)];
end

end