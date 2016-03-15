function [H, P, llr, exitflag] = circ_vmum_test_ts_k(data_x, data_y, alpha, options)
%CIRC_VMUM_TEST_SS one sample test for vMUM model
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
[mi_x_0, k_0, p1_x_0, p2_x_0, p3_x_0, ...
 mi_y_0, p1_y_0, p2_y_0, p3_y_0, ll_0, exitflag_0] = ...
    circ_vmum_test_ts_k_null(data_x, data_y, options);

%% Calculate ll of alternate hypothesis
[ll_1, exitflag_1] = circ_vmum_test_ts_k_alt(data_x, data_y, ...
                            mi_x_0, k_0, p1_x_0, p2_x_0, p3_x_0, ...
                            mi_y_0, p1_y_0, p2_y_0, p3_y_0, options);

%% Run test
llr = 2*(ll_1-ll_0);
P = 1 - chi2cdf(llr, 1);
H = P < alpha;
exitflag = min(exitflag_0, exitflag_1);

end

function [mi_x_0, k_0, p1_x_0, p2_x_0, p3_x_0, ...
          mi_y_0, p1_y_0, p2_y_0, p3_y_0, ll_0, exitflag_0] = ...
    circ_vmum_test_ts_k_null(data_x, data_y, options)

%% Calculate ML for the entire set
[mi_ml, k_ml, p1_ml, p2_ml, p3_ml] = ...
            circ_vmum_est_ml([data_x; data_y],...
                             [], [], [], [], [], ...
                             options);

%% Calculate ML for separate set using parameters of entire set as initial point
[params_a_0, ll_a_0_neg, exitflag_a_0] = fmincon(@(params) ...
                   -circ_vmum_ll(params(1), params(2), params(3), ...
                                 params(4), params(5), data_x, true)...
                   -circ_vmum_ll(params(6), params(2), params(7), ...
                                 params(8), params(9), data_y, true), ...
                   [mi_ml, k_ml, p1_ml, p2_ml, p3_ml...
                    mi_ml, p1_ml, p2_ml, p3_ml], [], [], ...
                   [0, 0, 1, 1, 1, 0, 0, 0, 0;...
                    0, 0, 0, 0, 0, 0, 1, 1, 1], [1; 1], ...
                   [-inf, -inf, 0, 0, 0, -inf, 0, 0, 0], ...
                   [inf, inf, 1, 1, 1, inf, 1, 1, 1], ...
                   [], options);
ll_a_0 = -ll_a_0_neg;      

%% Calculate ML for separate set using parameters of entire set as initial point
[mi_x_ml, ~, p1_x_ml, p2_x_ml, p3_x_ml] = circ_vmum_est_ml(data_x);
[mi_y_ml, ~, p1_y_ml, p2_y_ml, p3_y_ml] = circ_vmum_est_ml(data_y);

[params_b_0, ll_b_0_neg, exitflag_b_0] = fmincon(@(params) ...
                   -circ_vmum_ll(params(1), params(2), params(3), ...
                                 params(4), params(5), data_x, true)...
                   -circ_vmum_ll(params(6), params(2), params(7), ...
                                 params(8), params(9), data_y, true), ...
                   [mi_x_ml, k_ml, p1_x_ml, p2_x_ml, p3_x_ml...
                    mi_y_ml, p1_y_ml, p2_y_ml, p3_y_ml], [], [], ...
                   [0, 0, 1, 1, 1, 0, 0, 0, 0;...
                    0, 0, 0, 0, 0, 0, 1, 1, 1], [1; 1], ...
                   [-inf, -inf, 0, 0, 0, -inf, 0, 0, 0], ...
                   [inf, inf, 1, 1, 1, inf, 1, 1, 1], ...
                   [], options);
ll_b_0 = -ll_b_0_neg;

%% Calculate which one is better
if ll_a_0 > ll_b_0
    params_0 = params_a_0;
    ll_0 = ll_a_0;
    exitflag_0 = exitflag_a_0;
else
    params_0 = params_b_0;
    ll_0 = ll_b_0;
    exitflag_0 = exitflag_b_0;
end

mi_x_0 = params_0(1);
k_0 = params_0(2);
p1_x_0 = params_0(3);
p2_x_0 = params_0(4);
p3_x_0 = params_0(5);
mi_y_0 = params_0(6);
p1_y_0 = params_0(7);
p2_y_0 = params_0(8);
p3_y_0 = params_0(9);

end



function [ll_1, exitflag_1] = ...
    circ_vmum_test_ts_k_alt(data_x, data_y, ...
                            mi_x_0, k_0, p1_x_0, p2_x_0, p3_x_0, ...
                            mi_y_0, p1_y_0, p2_y_0, p3_y_0, options)
% Calculate ML estimate of parameters, using the previous ML
% estimate as starting point. This ensures that LLRT > 0, 
% but is not guaranteed to find the best log-likelihood, 
% since it often happens that algorithm gets stuck 
% in a local minimum near the solution with mu=mu_0.
[~, ~, ~, ~, ~, ll_x_a, exitflag_x_a, ~] = ...
    circ_vmum_est_ml(data_x, mi_x_0, k_0, p1_x_0, p2_x_0, p3_x_0, options);
[~, ~, ~, ~, ~, ll_y_a, exitflag_y_a, ~] = ...
    circ_vmum_est_ml(data_y, mi_y_0, k_0, p1_y_0, p2_y_0, p3_y_0, options);
ll_a = ll_x_a + ll_y_a;
exitflag_a = min(exitflag_x_a, exitflag_y_a);

% Calculate ML estimate of mu_hat, k_hat and p, using the vMUM-MM
% estimate as starting point. This usually finds higher log-likelihood.
% However, since this is not guaranteed, at the end we still compare 
% with the result of the previous one, and choose the highest.
[~, ~, ~, ~, ~, ll_x_b, exitflag_x_b] = ...
            circ_vmum_est_ml(data_x, [], [], [], [], [], options);
[~, ~, ~, ~, ~, ll_y_b, exitflag_y_b] = ...
            circ_vmum_est_ml(data_y, [], [], [], [], [], options);
ll_b = ll_x_b + ll_y_b;
exitflag_b = min(exitflag_x_b, exitflag_y_b);
        
% Choose the higher one
if ll_a > ll_b
    ll_1 = ll_a;
    exitflag_1 = exitflag_a;
else
    ll_1 = ll_b;
    exitflag_1 = exitflag_b;
end

end

