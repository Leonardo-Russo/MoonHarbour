function perturbation = get_state_perturbation(state)
% Description: the state is assumed to be position and velocity.

n = length(state);
% delta = 0.05;
delta = 0.1;
perturbation = zeros(n, 1);

for k = 1 : n
    mu = 0;
    sigma = abs(delta*state(k));
    truncation = 3*sigma;

    dist = makedist('normal', 'mu', mu, 'sigma', sigma);
    dist = truncate(dist, -truncation, truncation);
    perturbation(k) = random(dist, 1);
end


end