function [mu, k1, k2, p1, p2, p3] = circ_vmumfb_standard(mu, k1, k2, p1, p2, p3)
%CIRC_VMUM_STANDARD   converts parameters of vMUMFB to standard
%   [mu, k1, k2, p1, p2, p3] = CIRC_VMUMFB_STANDARD(mu, k1, k2, p1, p2, p3) takes 
%   the parameters of the vMUMFB distribution to a standard form 
%   with k1>0 and k2>0, p1 >= p2, and mu between 0 and 2*pi
% 
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

%% Asserts
circ_vmumfb_asserts(mu, k1, k2, p1, p2, p3);

%% Convert k<0
if k1 < 0 && k2 < 0
    k1 = -k1;
    k2 = -k2;
    mu = mu + pi;
end
if (k1 < 0 && k2 > 0) || (k1 > 0 && k2 < 0)
    warning('k1 and k2 have opposite sign, which means that this does not satisfy the front/back property of vMUMFB. ')
end

%% Convert p1 < p2
if p1 < p2
    mu = pi - mu;
    temp = p1;
    p1 = p2;
    p2 = temp;
    
    temp = k1;
    k1 = k2;
    k2 = temp;
end
    
%% Convert mu
mu = mod(mu, 2*pi);

%% Correct p1, p2, p3 (sometimes these numbers deviate slightly from sum=1)
temp = p1 + p2 + p3;
p1 = p1./temp;
p2 = p2./temp;
p3 = p3./temp;
