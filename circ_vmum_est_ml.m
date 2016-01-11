function [mu_hat, k_hat, p1_hat, p2_hat, p3_hat, ll, exitflag, output] = ...
    circ_vmum_est_ml(data, mu_0, k_0, p1_0, p2_0, p3_0, options)
%CIRC_VMUM_EST_ML    calculates ML estimate of vMUM model parameters
%   CIRC_VMUM_EST_ML(data) calculates ML estimate of vMUM model parameters
%   with initial point given by the vMUM-MME estimator
%
%   CIRC_VMUM_EST_ML(data, mu_0, k_0, p1_0, p2_0, p3_0) calculates
%   estimate using [mu_0, k_0, p1_0, p2_0, p3_0] as initial point
%
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena


%% Initial point
if nargin == 1
    [mu_0, k_0, p1_0, p2_0, p3_0] = circ_vmum_est_mm(data);
end

%% Add options to fmincon
if nargin < 7
    options = optimoptions('fmincon', ...
        'Display', 'none', ...
        'Algorithm', 'sqp');
end

%% Asserts
assert(p1_0>=0.0 & p1_0<=1.0);
assert(p2_0>=0.0 & p2_0<=1.0);
assert(p3_0>=0.0 & p3_0<=1.0);
assert(abs((p1_0+p2_0+p3_0)-1.0)<1E-10);
assert(isscalar(mu_0) & isscalar(k_0));
assert(isscalar(p1_0) & isscalar(p2_0) & isscalar(p3_0));
assert(isvector(data));


%% Run
[params, ll_neg, exitflag, output] = fmincon(@(params) ...
                   -circ_vmum_ll(params(1), params(2), params(3), ...
                                 params(4), params(5), data), ...
                   [mu_0, k_0, p1_0, p2_0, p3_0], [], [], ...
                   [0, 0, 1, 1, 1], 1, ...
                   [-inf, -inf, 0, 0, 0], [inf, inf, 1, 1, 1], ...
                   [], options);
       
if exitflag <= 0
   warning(strcat('exitflag was ', num2str(exitflag)));
end

%% Assign output values
ll = -ll_neg;

[mu_hat, k_hat, p1_hat, p2_hat, p3_hat] = ...
    circ_vmum_standard(params(1), params(2), ...
                       params(3), params(4), params(5));

end
