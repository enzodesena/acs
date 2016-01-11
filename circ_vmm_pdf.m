function pdf = circ_vmm_pdf(mu, k, p, theta)
%CIRC_VMM_PDF returns PDF values of von Mises mixture
% 
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

%% Asserts
assert(isvector(mu));
assert(isvector(k));
assert(isvector(p));
assert(length(mu) == length(k));
assert(length(mu) == length(p));
assert(isvector(theta));

pdf = zeros(size(theta));
for i=1:length(mu)
    pdf = pdf + p(i).*circ_vm_pdf(mu(i), k(i), theta);

end

assert(all(not(isnan(pdf))));
