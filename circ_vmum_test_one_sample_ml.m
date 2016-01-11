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
        'Display', 'none', ...
        'Algorithm', 'sqp');
[~, ll_0_neg, exitflag, ~] = fmincon(@(params) ...
                   -circ_vmum_ll(mu, k_mm_0, p1_mm_0, p2_mm_0, p3_mm_0, data), ...
                   [k_mm_0, p1_mm_0, p2_mm_0, p3_mm_0], [], [], ...
                   [0, 1, 1, 1], 1, ...
                   [-inf, 0, 0, 0], [inf, 1, 1, 1], ...
                   [], options);
       
if exitflag <= 0
   warning(strcat('exitflag was ', num2str(exitflag)));
end

ll_0 = -ll_0_neg;

%% Calculate ll of alternate hypothesis
[mu_mm_1, k_mm_1, p1_mm_1, p2_mm_1, p3_mm_1] = circ_vmum_est_mm(data);
[~, ~, ~, ~, ~, ll_1, ~, ~] = ...
    circ_vmum_est_ml(data, mu_mm_1, k_mm_1, p1_mm_1, p2_mm_1, p3_mm_1, options);
   

%% Run test
llrt = ll_1 - ll_0;
P = 1 - chi2cdf(2*llrt, 1);
H = P < alpha;

end
