function pdf = circ_vm_pdf(mu, k, theta, large_k_correction)
%CIRC_VM_PDF von Mises distribution PDF
%   pdf = CIRC_VM_PDF(mu, k, theta) returns the values of the PDF
%   of the von Mises distribution with parameters mu and k
%   at angle(s) theta
%   
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

if nargin == 3
    large_k_correction = true;
end

%% Asserts
assert(isscalar(mu));
assert(isscalar(k));
assert(isvector(theta));

if k<0
    k = -k;
    mu = mu + pi;
end

%% Run

pdf = 1./(2.*pi.*besseli(0,k)).*exp(k.*cos(theta-mu));

%% Large k correction
% For large k the above expression given nan. 
% Using normal approximation
if large_k_correction && abs(k) > 500
    pdf = normpdf(angledif(theta, mu).*sqrt(k)).*sqrt(k);
end

assert(all(not(isnan(pdf))));
