function [H, P, llrt] = circ_vm_test_ss(data, mu_0, alpha)
%CIRC_VMUM_REST_ONE_SAMPLE_MM one sample test for vMUM model
%
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena


if nargin <= 1 || isempty(mu_0)
    mu_0 = 0;
end

if nargin <= 2 || isempty(alpha)
    alpha = 0.05;
end


%% Assert
assert(iscolumn(data));
assert(isscalar(mu_0));
assert(isscalar(alpha));

%% Run test
N = length(data);
[~, k_hat] = circ_vm_est(data);
[~, k_hat_0] = circ_vm_est(data, mu_0);

data_circ_r = circ_r(data);

C_0 = mean(cos(data - mu_0));
llrt = 2*N*(k_hat*data_circ_r - k_hat_0*C_0 - log(besseli(0, k_hat)) + log(besseli(0, k_hat_0)));

P = 1 - chi2cdf(llrt, 1);
H = P < alpha;

end
