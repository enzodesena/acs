function circ_vmumfb_asserts(mu, k1, k2, p1, p2, p3, tollerance)
%CIRC_VMUM_ASSERTS checks consistency of vMUM parameters
%   
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

if nargin == 6
    tollerance = 1e-4;
end

assert(isscalar(mu) & isscalar(k1) & isscalar(k2));
assert(isscalar(p1) & isscalar(p2) & isscalar(p3));
assert(p1>-tollerance & p1<1.0+tollerance);
assert(p2>-tollerance & p2<1.0+tollerance);
assert(p3>-tollerance & p3<1.0+tollerance);
if abs((p1+p2+p3)-1.0)>tollerance
    display([p1, p2, p3, p1+p2+p3])
    assert(false);
end
assert(abs((p1+p2+p3)-1.0)<3*tollerance);

   