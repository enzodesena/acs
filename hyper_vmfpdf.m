function pdf = hyper_vmfpdf(mu, k, x)
% pdf = hyper_vmfpdf(mus, kappa, x)
% PDF value of a von-Mises-Fisher distribution
% mu: 1xD
% k:  1x1 scalar value of the concentraiton parameter
% x: 1xD or NxD, in which case the output is Nx1, containing
%    the relative pdf values in each row.


%% Asserts
assert(isvector(mu));
mu = mu(:);
assert(abs(mu'*mu-1)<10*eps);
[N, D] = size(x);
assert(length(mu) == D);
assert(isscalar(k) & k>=0);

%% Calculate
constant = (k^(D/2-1))/((2*pi)^(D/2)*besseli(D/2-1, k));
pdf = nan(N, 1);
for n=1:N
    pdf(n) = exp(k.*mu'*x(n, :)').*constant;
end
end