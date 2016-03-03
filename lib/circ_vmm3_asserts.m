function circ_vmm3_asserts(mu, k, p, tol)
%CIRC_VMM3_ASSERTS checks consistency of vMM3 parameters
%   
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

if nargin == 3 || isempty(tol)
    tol = 1e-4;
end

assert(isscalar(mu) & isscalar(k) & isscalar(p));
assert(p>-tol & p<1.0+tol);

   