function [mus, kappas, alphas, posterior] = ...
    hyper_vmfm_est_em(x, num_components, num_iterations, ...
                      alphas, mus, kappas)
%hyper_vmfm_est_em 
% 
%   Directional Statistics Toolbox (DirStat)
%   Copyright 2016 Enzo De Sena

[N, D] = size(x);
assert(D == 3);


if nargin <= 2
    num_iterations = 10;
end
if nargin <= 3
    alphas = ones(num_components,1)/num_components;
end
if nargin <= 4
    mus = rand(num_components, D);
    mus = mus ./ repmat(sqrt(sum(mus.^2, 2)), 1, D);
end
if nargin <= 5
    kappas = rand(num_components, 1);
end

%% Initialization
posterior = nan(N, num_components);

for j=1:num_iterations
    ptot = zeros(N, 1);
    for h=1:num_components
        ptot = ptot + alphas(h).*hyper_vmfpdf(mus(h, :), kappas(h), x);
    end

    % ps is a matrix that contains the probability 
    
    for h=1:num_components
        posterior(:, h) = alphas(h).*hyper_vmfpdf(mus(h, :), kappas(h), x)./ptot;
        
        alphas(h) = mean(posterior(:, h));
        rh = sum(repmat(posterior(:, h), 1, D).*x, 1);
        mus(h, :) = rh./norm(rh);
        kappas(h) = fsolve(@(kappa)besseli(D/2,kappa)./besseli(D/2-1,kappa)-norm(rh)/(N*alphas(h)), 1);
    end
end
end
