function [mu, k, p] = circ_vmm3_standard(mu, k, p)
%CIRC_VMM3_STANDARD   converts parameters of vMUM to standard
%   [mu, k, p] = CIRC_VMM3_STANDARD(mu, k, p) takes 
%   the parameters of the vMM3 distribution to a standard form 
%   with k>0, p >= 1-p, and mu between 0 and 2*pi
% 
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

%% Asserts
circ_vmm3_asserts(mu, k, p);

%% Convert k<0
if k < 0
    k = -k;
    mu = mu + pi;
end

%% Convert p < 1-p
if p < 1-p;
    mu = mu + pi;
    p = 1-p;
end
    
%% Convert mu
mu = mod(mu, 2*pi);
