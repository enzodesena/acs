function theta = circ_vmumfb_rnd(mu, k1, k2, p1, p2, p3, N)
%CIRC_VMUM_PDF returns PDF values of vMUM model
% 
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

%% Asserts
assert(isscalar(mu));
assert(isscalar(k1));
assert(isscalar(k2));
assert(isscalar(p1));
assert(isscalar(p2));
assert(isscalar(p3));
assert(isscalar(N));

%% Run
theta = circ_vmm_rnd([mu, pi-mu, 0], [k1, k2, 0], [p1, p2, p3], N);