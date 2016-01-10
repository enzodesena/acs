function [mu_hat, k_hat] = circ_vm_est(data, mu)
%CIRC_VM_EST maximum likelihood estimates of von Mises data
%   [mu_hat, k_hat] = CIRC_VM_EST(data) returns the ML estimate of mu_hat, 
%   k_hat for von Mises distribution.
% 
%   [mu_hat, k_hat] = CIRC_VM_EST(data, mu) returns 
%   the estimate of k_hat under the prior that mu_hat=mu
% 
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

data_circ_mean = circ_mean(data);
data_circ_r = circ_r(data);

if nargin < 2
    mu_hat = data_circ_mean;
else
    mu_hat = mu;
end

options = optimset('Display', 'off') ;
k_hat = fsolve(@(x)besseli(1,x)./besseli(0,x) - ...
        data_circ_r.*cos(data_circ_mean-mu_hat), 1, options);

