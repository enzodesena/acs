function [H, P, llrt] = circ_vm_test_ts_mu_d(data_x, data_y, alpha)
%CIRC_VMUM_REST_ONE_SAMPLE_MM two sample test for vMUM model
% with different estimates of 'k'. I derived this test, but I need to
% double-check the math.
%
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena



if nargin <= 2 || isempty(alpha)
    alpha = 0.05;
end


%% Assert
assert(iscolumn(data_x));
assert(iscolumn(data_y));
assert(isscalar(alpha));

%% Run test
N_x = length(data_x);
N_y = length(data_y);
N = N_x + N_y;
R_x = circ_r(data_x);
R_y = circ_r(data_y);
R = circ_r([data_x; data_y]);


theta_ol = circ_mean([data_x; data_y]);
theta_x_ol = circ_mean(data_x);
theta_y_ol = circ_mean(data_y);

mi_tilde = theta_ol;
[~, k_x_hat] = circ_vm_est(data_x);
[~, k_y_hat] = circ_vm_est(data_y);
[~, k_x_tilde] = circ_vm_est(data_x, mi_tilde);
[~, k_y_tilde] = circ_vm_est(data_y, mi_tilde);

llrt = 2*N_x*(k_x_hat*R_x-log(besseli(0, k_x_hat))) ...
    + 2*N_y*(k_y_hat*R_y-log(besseli(0, k_y_hat))) ...
    - 2*N_x*(k_x_tilde*R_x*cos(theta_x_ol-theta_ol)-log(besseli(0, k_x_tilde))) ...
    - 2*N_y*(k_y_tilde*R_y*cos(theta_y_ol-theta_ol)-log(besseli(0, k_y_tilde)));

P = 1 - chi2cdf(llrt, 1);
H = P < alpha;

end
