function theta = circ_vmum_rnd(mu, k, p1, p2, p3, N)
%CIRC_VMUM_PDF returns PDF values of vMUM model
% 
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

%% Asserts
assert(isscalar(mu));
assert(isscalar(k));
assert(isscalar(p1));
assert(isscalar(p2));
assert(isscalar(p3));
assert(isscalar(N));

%% Run
theta = circ_vmm_rnd([mu, mu+pi, 0], [k, k, 0], [p1, p2, p3], N);