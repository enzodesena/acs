function [H, P, llrt, exitflag] = circ_vmum_test_ss(data, mi, alpha, options)
%CIRC_VMUM_TEST_SS one sample test for vMUM model
%
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena


if nargin <= 1 || isempty(mi)
    mi = 0;
end

if nargin <= 2 || isempty(alpha)
    alpha = 0.05;
end

if nargin <= 3    
    options = optimoptions('fmincon', ...
            'Display', 'notify-detailed', ...
            'Algorithm', 'sqp', ...
            'MaxFunEvals', 2000);
end

%% Assert
assert(iscolumn(data));
assert(isscalar(mi));
assert(isscalar(alpha));

%% Calculate ll of null hypothesis
% Calculate starting point using vMUM-MM method, assuming mu=mu_0
[mu_mm_0, k_mm_0, p1_mm_0, p2_mm_0, p3_mm_0] = circ_vmum_est_mm(data, mi);
assert(abs(mu_mm_0-mi)<1e-10);

% Calculate ML of parameters, assuming mu=mu_0
[params_0, ll_0_neg, exitflag_0] = fmincon(@(params) ...
                   -circ_vmum_ll(mi, params(1), params(2), ...
                                 params(3), params(4), data, true), ...
                   [k_mm_0, p1_mm_0, p2_mm_0, p3_mm_0], [], [], ...
                   [0, 1, 1, 1], 1, ...
                   [-inf, 0, 0, 0], [inf, 1, 1, 1], ...
                   [], options);
ll_0 = -ll_0_neg;

%% Calculate ll of alternate hypothesis
% Calculate ML estimate of parameters, using the previous ML
% estimate as starting point. This ensures that LLRT > 0, 
% but is not guaranteed to find the best log-likelihood, 
% since it often happens that algorithm gets stuck 
% in a local minimum near the solution with mu=mu_0.
[~, ~, ~, ~, ~, ll_1_a, exitflag_1_a, ~] = ...
    circ_vmum_est_ml(data, mi, params_0(1), ...
                     params_0(2), params_0(3), params_0(4), options);
                 
% Calculate ML estimate of mu_hat, k_hat and p, using the vMUM-ML
% estimate as starting point. This usually finds higher log-likelihood.
% However, since this is not guaranteed, at the end we still compare 
% with the result of the previous one, and choose the highest.
[mu_mm, k_mm, p1_mm, p2_mm, p3_mm] = circ_vmum_est_mm(data);
[~, ~, ~, ~, ~, ll_1_b, exitflag_1_b, ~] = ...
    circ_vmum_est_ml(data, mu_mm, k_mm, p1_mm, p2_mm, p3_mm, options);

% Choose the higher one
if ll_1_a > ll_1_b
    ll_1 = ll_1_a;
    exitflag_1 = exitflag_1_a;
else
    ll_1 = ll_1_b;
    exitflag_1 = exitflag_1_b;
end

%% Run test
llrt = 2*(ll_1-ll_0);
P = 1 - chi2cdf(llrt, 1);
H = P < alpha;
exitflag = min(exitflag_0, exitflag_1);

end
