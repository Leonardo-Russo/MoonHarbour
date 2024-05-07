function [beta_dist, gamma_dist] = define_distributions(MC)

% Define gamma as a normally distributed variable
sigma_gamma = deg2rad(2);           % gamma standard deviation
mu_gamma = 0;                       % gamma mean value
truncation_gamma = 3*sigma_gamma;   % gamma truncation limit

gamma_dist = makedist('normal', 'mu', mu_gamma, 'sigma', sigma_gamma);
gamma_dist = truncate(gamma_dist, -truncation_gamma, truncation_gamma);
gamma_dist = random(gamma_dist, MC, 1);


% Generate beta as a uniform random variable between 0 and 2*pi
beta_lower = 0;
beta_upper = 2*pi;
beta_dist = makedist('uniform', 'Lower', beta_lower, 'Upper', beta_upper);
beta_dist = random(beta_dist, MC, 1);

end