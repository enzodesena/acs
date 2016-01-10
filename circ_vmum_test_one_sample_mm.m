function [H, P, llrt] = circ_vmum_test_one_sample_mm(data, mu)
%CIRC_VMUM_REST_ONE_SAMPLE_MM one sample test for vMUM model
%
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena


alpha = 0.05;

if nargin == 1
    mu = 0;
end

%% Assert
assert(isvector(data));
assert(isscalar(mu));


[mu_hat_0, k_hat_0, p1_hat_0, p2_hat_0, p3_hat_0] = circ_vmum_est_mm(data, mu);
[mu_hat_1, k_hat_1, p1_hat_1, p2_hat_1, p3_hat_1] = circ_vmum_est_mm(data);

llrt = getll(mu_hat_1, k_hat_1, p1_hat_1, p2_hat_1, p3_hat_1, data) - ...
       getll(mu_hat_0, k_hat_0, p1_hat_0, p2_hat_0, p3_hat_0, data);

P = 1 - chi2cdf(2*llrt, 1);

H = P < alpha;

end

function l = getll(mu, kappa, p1, p2, p3, data)
l = sum(log(circ_vmm_pdf([mu, mu+pi, 0], ...
                         [kappa, kappa, 0], ...
                         [p1, p2, p3], data)));
end
