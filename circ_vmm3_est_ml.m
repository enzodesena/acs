function [mu_hat, k_hat, p_hat, ll, exitflag, output] = ...
    circ_vmm3_est_ml(data, mu_0, k_0, p_0, options)
%CIRC_VMM3_EST_ML    calculates ML estimate of vMM3 model parameters
%   CIRC_VMM3_EST_ML(data) calculates ML estimate of vMM3 model parameters
%   with initial point given by the vMUM-MME estimator
%
%   CIRC_VMM3_EST_ML(data, mu_0, k_0, p_0) calculates
%   estimate using [mu_0, k_0, p_0] as initial point
%
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena


%% Initial point
assert(iscolumn(data));
if nargin == 1
    [mu_0, k_0, p_0] = circ_vmm3_est_mm(data);
end

%% Add options to fmincon
if nargin < 7
    options = optimoptions('fmincon', ...
        'Display', 'notify-detailed', ...
        'Algorithm', 'sqp');
end

%% Asserts
assert(isscalar(mu_0) & isscalar(k_0) & isscalar(p_0));
assert(isvector(data));

%% Run
[mu_0, k_0, p_0] = circ_vmm3_standard(mu_0, k_0, p_0);
[params, ll_neg, exitflag, output] = fmincon(@(params) ...
                   -circ_vmm3_ll(mu_0, k_0, p_0, data), ...
                   [mu_0, k_0, p_0], [], [], ...
                   [], [], ...
                   [-inf, -inf, 0], [inf, inf, 1], ...
                   [], options);
       

%% Assign output values
ll = -ll_neg;

[mu_hat, k_hat, p_hat] = circ_vmm3_standard(params(1), params(2), params(3));

end
