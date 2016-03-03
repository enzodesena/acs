function [H, P, llrt, exitflag] = circ_vmum_test_ts(data_a, data_b, alpha, options)
%CIRC_VMUM_TEST_SS one sample test for vMUM model
%
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

assert(iscolumn(data_a));
assert(iscolumn(data_b));

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
assert(iscolumn(data_a));
assert(iscolumn(data_b));
assert(isscalar(alpha));

%% Calculate ll of null hypothesis
% Calculate mu for entire set
mu = circ_mean([data_a(:); data_b(:)]);

% Calculate ML of parameters, assuming mu=mu_0 and with 
% starting point using vMUM-MM method
[k_a_ml_0, p1_a_ml_0, p2_a_ml_0, p3_a_ml_0, ll_a_0, exitflag_a_0] = ...
    circ_vmum_est_ml_null(data_a, mu, options);
[k_b_ml_0, p1_b_ml_0, p2_b_ml_0, p3_b_ml_0, ll_b_0, exitflag_b_0] = ...
    circ_vmum_est_ml_null(data_b, mu, options);
ll_0 = ll_a_0 + ll_b_0;
exitflag_0 = min(exitflag_a_0, exitflag_b_0);

%% Calculate ll of alternate hypothesis
[ll_a_1, exitflag_a_1] = ...
    circ_vmum_est_ml_alt(data_a, mu, k_a_ml_0, p1_a_ml_0, p2_a_ml_0, p3_a_ml_0, options);
[ll_b_1, exitflag_b_1] = ...
    circ_vmum_est_ml_alt(data_b, mu, k_b_ml_0, p1_b_ml_0, p2_b_ml_0, p3_b_ml_0, options);
ll_1 = ll_a_1 + ll_b_1;
exitflag_1 = min(exitflag_a_1, exitflag_b_1);

%% Run test
llrt = 2*(ll_1-ll_0);
P = 1 - chi2cdf(llrt, 1);
H = P < alpha;
exitflag = min(exitflag_0, exitflag_1);

end

function [k_ml_0, p1_ml_0, p2_ml_0, p3_ml_0, ll_0, exitflag_0] = ...
    circ_vmum_est_ml_null(data, mu, options)
[mu_mm_0, k_mm_0, p1_mm_0, p2_mm_0, p3_mm_0] = ...
    circ_vmum_est_mm(data, mu);
assert(abs(mu_mm_0-mu)<1e-10);
[params_ml_0, ll_0_neg, exitflag_0] = fmincon(@(params) ...
                   -circ_vmum_ll(mu, params(1), params(2), ...
                                 params(3), params(4), data, true), ...
                   [k_mm_0, p1_mm_0, p2_mm_0, p3_mm_0], [], [], ...
                   [0, 1, 1, 1], 1, ...
                   [-inf, 0, 0, 0], [inf, 1, 1, 1], ...
                   [], options);
               
k_ml_0 = params_ml_0(1);
p1_ml_0 = params_ml_0(2);
p2_ml_0 = params_ml_0(3);
p3_ml_0 = params_ml_0(4);
               
ll_0 = -ll_0_neg;

end



function [ll_1, exitflag_1] = ...
    circ_vmum_est_ml_alt(data, mu, k_ml_0, p1_ml_0, p2_ml_0, p3_ml_0, options)
% Calculate ML estimate of parameters, using the previous ML
% estimate as starting point. This ensures that LLRT > 0, 
% but is not guaranteed to find the best log-likelihood, 
% since it often happens that algorithm gets stuck 
% in a local minimum near the solution with mu=mu_0.
[~, ~, ~, ~, ~, ll_1_a, exitflag_1_a, ~] = ...
    circ_vmum_est_ml(data, mu, k_ml_0, p1_ml_0, p2_ml_0, p3_ml_0, options);

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
end