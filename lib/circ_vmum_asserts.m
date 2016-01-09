function circ_vmum_asserts(mu, k, p1, p2, p3)
%CIRC_VMUM_ASSERTS checks consistency of vMUM parameters
%   
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

assert(p1>=0.0 & p1<=1.0);
assert(p2>=0.0 & p2<=1.0);
assert(p3>=0.0 & p3<=1.0);
assert(abs((p1+p2+p3)-1.0)<1E-3);
assert(isscalar(mu) & isscalar(k));
assert(isscalar(p1) & isscalar(p2) & isscalar(p3));