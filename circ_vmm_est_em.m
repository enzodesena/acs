function circ_vmm_est_em(data)
%CIRC_VMM_EST_EM expectation maximisation estimation of von Mises mixtures
% 
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

data = data(:);

all_data = linspace(-pi, pi, 1000);

%% Calculate intial guess using MvM model
[mu_hat_MvM, k_hat_MvM, p_MvM] = circ_vmm3est(data);

%% Run EM

mu_hat = [mu_hat_MvM, mu_hat_MvM + pi, 0];
k_hat = [k_hat_MvM, k_hat_MvM, 0];
p_hat = [p_MvM, 1-p_MvM, 0.0001];

L = length(p_hat);

for i=1:100
    hold off
    circ_hist(data, 5/180*pi);
    
    hold on
    plot(all_data, circ_vmmpdf(mu_hat, k_hat, p_hat, all_data), 'LineWidth', 4);
    hold off
    pause(0.1)
    
    for l=1:L
        % Expectation
        p_mix = zeros(size(data));
        for m=1:L
            p_mix = p_mix + p_hat(m)*circ_vmpdf(mu_hat(m), k_hat(m), data);
        end
        p_cond = p_hat(l)* circ_vmpdf(mu_hat(l), k_hat(l), data)./p_mix;
        p_hat(l) = mean(p_cond);

        % Maximization
        r_l = [sum(cos(data).*p_cond), sum(sin(data).*p_cond)];
        mu_l = r_l./norm(r_l);
        mu_hat(l) = angle(mu_l(1) + 1i*mu_l(2));
        x_m = norm(r_l)/sum(p_cond);
        k_hat(l) = fsolve(@(x)besseli(1,x)./besseli(0,x)-x_m, 1);
    end

end
end

