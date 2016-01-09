function ll = circ_vmum_ll(mu, k, p1, p2, p3, data)
%CIRC_VMUM_LL returs log-likelihood of the vMUM model
%
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

%% Asserts
circ_vmum_asserts(mu, k, p1, p2, p3);

%% Calculate
ll = sum(log(circ_vmm_pdf([mu, mu+pi, 0], [k, k, 0], [p1, p2, p3], data)));

end
