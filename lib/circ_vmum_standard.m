function [mu, k, p1, p2, p3] = circ_vmum_standard(mu, k, p1, p2, p3)
%CIRC_VMUM_STANDARD   converts parameters of vMUM to standard
%   [mu, k, p1, p2, p3] = CIRC_VMUM_STANDARD(mu, k, p1, p2, p3) takes 
%   the parameters of the vMUM distribution to a standard form 
%   with k>0, p1 >= p2, and mu between 0 and 2*pi
% 
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

%% Asserts
circ_vmum_asserts(mu, k, p1, p2, p3);

%% Convert k<0
if k < 0
    k = -k;
    mu = mu + pi;
end

%% Convert p1 < p2
if p1 < p2
    mu = mu + pi;
    temp = p1;
    p1 = p2;
    p2 = temp;
end
    
%% Convert mu
mu = mod(mu, 2*pi);
