function [H, P, llr, exitflag] = circ_vmum_test_ts_mu(data_x, data_y, alpha, options)
%CIRC_VMUM_TEST_SS two sample test for vMUM model
%   assuming they have the same concentration parameter.
%
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

assert(iscolumn(data_x));
assert(iscolumn(data_y));

if nargin <= 3 || isempty(alpha)
    alpha = 0.05;
end
assert(isscalar(alpha));

if nargin <= 4    
    options = optimoptions('fmincon', ...
            'Display', 'notify-detailed', ...
            'Algorithm', 'sqp', ...
            'MaxFunEvals', 2000);
end

%% Assert
assert(iscolumn(data_x));
assert(iscolumn(data_y));
assert(isscalar(alpha));

%% Calculate ll of null hypothesis
% Calculate ML of parameters, assuming mu_a=mu_b and with 
% starting point using vMUM-MM method
[mu_ml_0, k_ml_0, p1_ml_0, p2_ml_0, p3_ml_0, ll_0, exitflag_0] = ...
    circ_vmum_est_ml([data_x; data_y]);

%% Calculate ll of alternate hypothesis
[~, ll_a_1_neg, exitflag_a_1] = fmincon(@(params) ...
                   -circ_vmum_ll(params(1), params(2), params(3), ...
                                 params(4), params(5), data_x, true) ...
                   -circ_vmum_ll(params(6), params(2), params(3), ...
                                 params(4), params(5), data_y, true), ...
                   [mu_ml_0, k_ml_0, p1_ml_0, p2_ml_0, p3_ml_0, mu_ml_0], ...
                   [], [], ...
                   [0, 0, 1, 1, 1, 0], 1, ...
                   [-inf, -inf, 0, 0, 0, -inf], [inf, inf, 1, 1, 1, inf], ...
                   [], options);
ll_a_1 = -ll_a_1_neg;

[~, ll_b_1_neg, exitflag_b_1] = fmincon(@(params) ...
                   -circ_vmum_ll(params(1), params(2), params(3), ...
                                 params(4), params(5), data_x, true) ...
                   -circ_vmum_ll(params(6), params(2), params(3), ...
                                 params(4), params(5), data_y, true), ...
                   [circ_mean(mod(data_x*2, 2*pi))/2, ...
                    k_ml_0, p1_ml_0, p2_ml_0, p3_ml_0, ...
                    circ_mean(mod(data_y*2, 2*pi))/2], ...
                   [], [], ...
                   [0, 0, 1, 1, 1, 0], 1, ...
                   [-inf, -inf, 0, 0, 0, -inf], [inf, inf, 1, 1, 1, inf], ...
                   [], options);
ll_b_1 = -ll_b_1_neg;

if ll_a_1 > ll_b_1
    ll_1 = ll_a_1;
    exitflag_1 = exitflag_a_1;
else
    ll_1 = ll_b_1;
    exitflag_1 = exitflag_b_1;
end

%% Run test
llr = 2*(ll_1-ll_0);
P = 1 - chi2cdf(llr, 1);
H = P < alpha;
exitflag = min(exitflag_0, exitflag_1);

end
