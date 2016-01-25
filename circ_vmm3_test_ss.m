function [H, P, llrt, exitflag] = circ_vmm3_test_ss(data, mu_0, alpha, options)
%CIRC_VMM3_TEST_SS one sample test for vMUM model
%
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena


if nargin <= 1 || isempty(mu_0)
    mu_0 = 0;
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
assert(isvector(data));
assert(isscalar(mu_0));
assert(isscalar(alpha));

%% Calculate ll of null hypothesis
% Calculate starting point using vMM3-MM method, assuming mu=mu_0
[mu_mm_0, k_mm_0, p_mm_0] = circ_vmm3_est_mm(data, mu_0);
assert(abs(mu_mm_0-mu_0)<1e-10);

% Calculate ML of k_hat and p, assuming mu=mu_0
[params_0, ll_0_neg, exitflag_0] = fmincon(@(params) ...
                   -sum(log(circ_vmm_pdf([mu_0, mu_0+pi], ...
                                         [params(1), params(1)], ...
                                         [params(2), 1-params(2)], data))), ...
                   [k_mm_0, p_mm_0], [], [], ...
                   [], [], ...
                   [-inf, 0], [inf, 1], ...
                   [], options);
ll_0 = -ll_0_neg;

%% Calculate ll of alternate hypothesis
% Calculate ML estimate of mu_hat, k_hat and p, using the previous ML
% estimate as starting poitn. This ensures that LLRT > 0, 
% but is not guaranteed to find the best log-likelihood, 
% since it often happens that algorithm gets stuck 
% in a local minimum near the solution with mu=mu_0.
[~, ll_1_neg_a, exitflag_1_a] = fmincon(@(params) ...
                   -sum(log(circ_vmm_pdf([params(1), params(1)+pi], ...
                                         [params(2), params(2)], ...
                                         [params(3), 1-params(3)], data))), ...
                   [mu_0, params_0(1), params_0(2)], [], [], ...
                   [], [], ...
                   [-inf, -inf, 0], [inf, inf, 1], ...
                   [], options);
ll_1_a = -ll_1_neg_a;
               
% Calculate ML estimate of mu_hat, k_hat and p, using the vMM3-ML
% estimate as starting point. This usually finds higher log-likelihood.
% However, since this is not guaranteed, at the end we still compare 
% with the result of the previous one, and choose the highest.
[mu_mm_1, k_mm_1, p_mm_1] = circ_vmm3_est_mm(data);
[~, ll_1_neg_b, exitflag_1_b] = fmincon(@(params) ...
   -sum(log(circ_vmm_pdf([params(1), params(1)+pi], ...
                         [params(2), params(2)], ...
                         [params(3), 1-params(3)], data))), ...
   [mu_mm_1, k_mm_1, p_mm_1], [], [], ...
   [], [], ...
   [-inf, -inf, 0], [inf, inf, 1], ...
   [], options);
ll_1_b = -ll_1_neg_b;

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
