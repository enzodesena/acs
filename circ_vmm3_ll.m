function ll = circ_vmm3_ll(mu, k, p, data)
%CIRC_VMUM_LL returs log-likelihood of the vMUM model
%
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

%% Asserts
assert(isscalar(mu));
assert(isscalar(k));
assert(isscalar(p));
assert(isvector(data));
circ_vmm3_asserts(mu, k, p);


%% Calculate
ll = sum(log(circ_vmm_pdf([mu, mu+pi], [k, k], [p, 1-p], data)));
assert(not(isnan(ll)));

end
