function ll = circ_vmum_ll(mu, k, p1, p2, p3, data)
%CIRC_VMUM_LL returs log-likelihood of the vMUM model
%
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

%% Asserts
assert(isscalar(mu));
assert(isscalar(k));
assert(isscalar(p1));
assert(isscalar(p2));
assert(isscalar(p3));
assert(isvector(data));
circ_vmum_asserts(mu, k, p1, p2, p3, 1e-5);

%% Calculate
ll = sum(log(circ_vmum_pdf(mu, k, p1, p2, p3, data)));
assert(not(isnan(ll)));

end
