function [mu_hat, k1_hat, k2_hat, p1_hat, p2_hat, p3_hat, ll, exitflag, output] = ...
    circ_vmumfb_est_ml(data, mu_0, k1_0, k2_0, p1_0, p2_0, p3_0, options)
%CIRC_VMUMFB_EST_ML    calculates ML estimate of vMUMFB model parameters
%   CIRC_VMUMFB_EST_ML(data) calculates ML estimate of vMUMFB model parameters
%   with initial point given by the vMUM-MME estimator
%
%   CIRC_VMUMFB_EST_ML(data, mu_0, k1_0, k2_0, p1_0, p2_0, p3_0) calculates
%   estimate using [mu_0, k1_0, k2_0, p1_0, p2_0, p3_0] as initial point
%
%   Audio Circular Statistics (ACS) library
%   Copyright 2020 Enzo De Sena


%% Initial point
assert(iscolumn(data));
if nargin <= 1 || isempty(mu_0) || isempty(k1_0) || isempty(k2_0) || ...
        isempty(p1_0) || isempty(p2_0) || isempty(p3_0)
    [mu_0_fs, k1_0_fs, k2_0_fs, p1_0_fs, p2_0_fs, p3_0_fs] = ...
        circ_vmumfb_est_fs(data, linspace(0, pi, 10), ...
                                 linspace(0, 100, 10), linspace(0, 100, 10), ...
                                 linspace(0, 1, 10), linspace(0, 1, 10));
    if nargin <= 1 || isempty(mu_0)
        mu_0 = mu_0_fs;
    end
    if nargin <= 2 || isempty(k1_0)
        k1_0 = k1_0_fs;
    end
    if nargin <= 3 || isempty(k2_0)
        k2_0 = k2_0_fs;
    end
    if nargin <= 4 || isempty(p1_0)
        p1_0 = p1_0_fs;
    end
    if nargin <= 5 || isempty(p2_0)
        p2_0 = p2_0_fs;
    end
    if nargin <= 6 || isempty(p3_0)
        p3_0 = p3_0_fs;
    end
end

%% Add options to fmincon
if nargin < 7
    options = optimoptions('fmincon', ...
        'Display', 'notify-detailed', ...
        'Algorithm', 'sqp');
end

%% Asserts
assert(isscalar(mu_0) & isscalar(k1_0) & isscalar(k2_0));
assert(isscalar(p1_0) & isscalar(p2_0) & isscalar(p3_0));
assert(isvector(data));
circ_vmumfb_asserts(mu_0, k1_0, k2_0, p1_0, p2_0, p3_0, 1e-6);


%% Run
[mu_0, k1_0, k2_0, p1_0, p2_0, p3_0] = ...
    circ_vmumfb_standard(mu_0, k1_0, k2_0, p1_0, p2_0, p3_0);
[params, ll_neg, exitflag, output] = fmincon(@(params) ...
                   -circ_vmumfb_ll(params(1), params(2), params(3), ...
                                   params(4), params(5), params(6), data, true, true), ...
                   [mu_0, k1_0, k2_0, p1_0, p2_0, p3_0], [], [], ...
                   [0, 0, 0, 1, 1, 1], 1, ...
                   [-inf, 0, 0, 0, 0, 0], [inf, 200, 200, 1, 1, 1], ...
                   [], options);
       

%% Assign output values
ll = -ll_neg;

[mu_hat, k1_hat, k2_hat, p1_hat, p2_hat, p3_hat] = ...
    circ_vmumfb_standard(params(1), params(2), params(3),...
                         params(4), params(5), params(6));
end
