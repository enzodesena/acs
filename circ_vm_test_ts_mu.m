function [H, P, llrt] = circ_vm_test_ts_mu(data_x, data_y, alpha)
%CIRC_VMUM_REST_ONE_SAMPLE_MM one sample test for vMUM model
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

[~, k_hat] = circ_vm_est([data_x; data_y]);

options = optimset('Display', 'off') ;
k_hat_12 = fsolve(@(x)besseli(1,x)./besseli(0,x) - ...
            (R_x.*N_x./N + R_y.*N_y./N), 1, options);
        
% Book: not working
%llrt = 2*(k_hat_12*(R_x+R_y) - k_hat*R - log(besseli(0, k_hat_12)) + log(besseli(0, k_hat)));
% Corrected from book, but only valid for N_x=N_y
%llrt = 2*N*(k_hat_12*(R_a+R_b)/2 - k_hat*R - log(besseli(0, k_hat_12)) + log(besseli(0, k_hat)));
llrt = 2*(k_hat_12*(N_x*R_x+N_y*R_y) - k_hat*N*R - N*log(besseli(0, k_hat_12)) + N*log(besseli(0, k_hat)));

P = 1 - chi2cdf(llrt, 1);
H = P < alpha;

end
