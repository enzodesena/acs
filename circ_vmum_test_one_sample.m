function [H, P, llrt] = circ_vmum_test_one_sample_ml(data, mu)
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

%% Calculate ll of null hypothesis
[~, k_mm_0, p1_mm_0, p2_mm_0, p3_mm_0] = circ_vmum_est_mm(data, mu);

options = optimoptions('fmincon', ...
        'Display', 'notify-detailed', ...
        'Algorithm', 'sqp');

[params, ll_0_neg] = fmincon(@(params) ...
                   -circ_vmum_ll(mu, params(1), params(2), params(3), params(4), data, true), ...
                   [k_mm_0, p1_mm_0, p2_mm_0, p3_mm_0], [], [], ...
                   [0, 1, 1, 1], 1, ...
                   [-inf, 0, 0, 0], [inf, 1, 1, 1], ...
                   [], options);
       
ll_0 = -ll_0_neg;

%% Calculate ll of alternate hypothesis
[~, ~, ~, ~, ~, ll_1, ~, ~] = ...
    circ_vmum_est_ml(data, mu, params(1), params(2), params(3), params(4), options);
   

%% Run test
llrt = -2*(ll_0-ll_1);
P = 1 - chi2cdf(2*llrt, 1);
H = P < alpha;

end
