function [mus, kappas, alphas, ps] = hyper_vmfm_est_em(x, H, J)
%CIRC_VMUM_TEST_SS two sample test for vMUM model
%   assuming they have the same concentration parameter.
%
%   Directional Statistics Toolbox (DirStat)
%   Copyright 2016 Enzo De Sena


[D, N] = size(x);

%% Initialization
alphas = ones(H,1)/H;
mus = rand(D, H);
mus = mus ./ repmat(sqrt(sum(mus.^2)), D, 1);
%mus = repmat([1; zeros(D-1,1)], 1, H);
kappas = rand(H, 1);

ps = nan(H, N);

for j=1:J

    ptot = zeros(N, 1);
    for h=1:H
        ptot = ptot + alphas(h)*sph_vmfpdf(mus(:, h), kappas(h), x);
    end

    % ps is a matrix that contains the probability 
    
    for h=1:H
        ps(h, :) = alphas(h)'.*sph_vmfpdf(mus(:, h), kappas(h), x)./ptot;
        
        alphas(h) = mean(ps(h, :));
        rh = sum(repmat(ps(h, :), D, 1).*x, 2);
        mus(:, h) = rh/norm(rh);
        kappas(h) = fsolve(@(x)besseli(D/2,x)./besseli(D/2-1,x)-norm(rh)/N, 1);
    end
    
    
    
end
    
end