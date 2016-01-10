function obs = circ_vmm_rnd(mu, kappa, p, N)
%CIRC_VMM_RND generates a random sample of a von Mises mixture
% 
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

%% Checks
assert(isvector(mu));
assert(isvector(kappa));
assert(isvector(p));
assert(length(mu) == length(kappa));
assert(length(mu) == length(p));

%% Definitions
obs = [];
L = length(mu);

%% Run
z_obs = drnd(1:L, p, N); % Observations of the hidden variable
for l=1:L
    obs = [obs; circ_vmrnd(mu(l), kappa(l), sum(z_obs == l))];
end
assert(length(obs) == N);

% Shuffle the elements
obs = obs(randperm(length(obs)));

end
