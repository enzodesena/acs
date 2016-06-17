function circ_vmum_asserts(mu, k, p1, p2, p3, tol)
%CIRC_VMUM_ASSERTS checks consistency of vMUM parameters
%   
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

if nargin == 5
    tol = 1e-4;
end

assert(isscalar(mu) & isscalar(k));
assert(isscalar(p1) & isscalar(p2) & isscalar(p3));
assert(p1>-tol & p1<1.0+tol);
assert(p2>-tol & p2<1.0+tol);
assert(p3>-tol & p3<1.0+tol);
if abs((p1+p2+p3)-1.0)>tol
    display([p1, p2, p3, p1+p2+p3])
    assert(false);
end
assert(abs((p1+p2+p3)-1.0)<3*tol);

   