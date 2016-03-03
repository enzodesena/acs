function pdf = circ_vmm3_pdf(mu_hat, k_hat, p, theta)
%CIRC_VMM3_PDF returns PDF values of the 3 parameters von Mises mixture
% 
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

pdf = 1./(2.*pi.*besseli(0,k_hat)).*(p.*exp(k_hat.*cos(theta-mu_hat))+...
    (1-p).*exp(-k_hat*cos(theta-mu_hat)));
