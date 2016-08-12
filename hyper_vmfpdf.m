function pdf = hyper_vmfpdf(mus, kappa, x)
%% Asserts
assert(isvector(mus));
mus = mus(:);
assert(abs(mus'*mus-1)<10*eps);
[D, N] = size(x);
assert(length(mus) == D);
assert(isscalar(kappa) & kappa>=0);

%% Calculate
c = (kappa^(D/2-1))/((2*pi)^(D/2)*besseli(D/2-1, kappa));
pdf = nan(N, 1);
for n=1:N
    pdf(n) = exp(kappa.*mus'*x(:, n)).*c;
end
end