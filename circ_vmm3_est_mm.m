function [mu_hat, k_hat, p] = circ_vmm3_est_mm(data, mu)
%CIRC_VMM3_EST_MM method of moment estimate of the vMM3 model
%   [mu_hat, k_hat, p] = CIRC_VMM3_EST_MM(data) returns the method 
%   of moment estimate of the parameters of the
%   three paramters of a mixture of von Mises distribution:
%   p*exp(k*cos(theta-mu)/I0(k) + (1-p)*exp(-k*cos(theta-mu)/I0(k)
%   using method explained in [Mardia and Jupp, 2000]
% 
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

%% Asserts
assert(iscolumn(data));

%% Estimation
psi = mod(data*2, 2*pi);

if nargin <= 1
    mu_hat = circ_mean(psi)/2;
else
    mu_hat = mu;
end

R2 = circ_r(psi);
options = optimset('Display', 'off') ;
k_hat = fsolve(@(x)1 - 2*besseli(1,x)./besseli(0,x)/x-R2, 1, options);

C = mean(cos(data));
S = mean(sin(data));
p = (besseli(1,k_hat)./besseli(0,k_hat) + C*cos(mu_hat) + S*sin(mu_hat))/2;
p = max(0, p);
p = min(1, p);

if k_hat<0
    k_hat = -k_hat;
    mu_hat = mu_hat + pi;
end
mu_hat = mod(mu_hat, 2*pi);
