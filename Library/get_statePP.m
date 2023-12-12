function PP = get_statePP(times, state)
% Description: this function retrieves the pp for each component of the
% state matrix intended as a MxN matrix, where N is the number of state
% components and M is the number of epochs.

N = size(state, 2);
PP = [];

for i = 1 : N
    
    PP = [PP; spline(times, state(:, i))];

end


end