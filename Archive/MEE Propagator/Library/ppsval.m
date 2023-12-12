function state = ppsval(PPs, time)
% Description: this function takes an array of pp's and computes the
% interpolated value for each component in the state.

N = length(PPs);
state = zeros(N, 1);

for i = 1 : N
    state(i) = ppval(PPs(i), time);
end


end