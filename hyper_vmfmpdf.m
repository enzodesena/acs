function pdf = hyper_vmfmpdf(mu, k, p, x)
%hyper_vmfmpdf returns the pdf values of 
% PDF value of a von-Mises-Fisher mixture distribution
% mu: HxD with H num of components in mixture and D dimension of the
% distributio
% k:  Hx1 values of the concentraiton parameter for each of the H component
% x: 1xD or NxD, in which case the output is Nx1, containing
%    the relative pdf values in each row.

%% Asserts
[H, D] = size(mu);
assert(isvector(k));
assert(length(k)==H);
k = k(:);
assert(isvector(p));
assert(length(p)==H);
p = p(:);
[N, Dx] = size(x);
assert(D==Dx);

%% Run
pdf = zeros(N, 1);

for h=1:H
    pdf = pdf + p(h).*hyper_vmfpdf(mu(h, :), k(h), x);
end

