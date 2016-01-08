function [H, P, llrt] = circ_vmum_test_one_sample(data, mu)

alpha = 0.05;

if nargin == 1
    mu = 0;
end

[mu_hat_0, k_hat_0, p1_hat_0, p2_hat_0, p3_hat_0] = circ_vmum_est_mm(data, mu);
[mu_hat_1, k_hat_1, p1_hat_1, p2_hat_1, p3_hat_1] = circ_vmum_est_mm(data);

llrt = getll(mu_hat_1, k_hat_1, p1_hat_1, p2_hat_1, p3_hat_1, data) - ...
       getll(mu_hat_0, k_hat_0, p1_hat_0, p2_hat_0, p3_hat_0, data);

P = 1 - chi2cdf(2*llrt, 1);

H = P < alpha;

end

function l = getll(mu, kappa, p1, p2, p3, data)
l = sum(log(circ_vmm_pdf([mu, mu+pi, 0], ...
                         [kappa, kappa, 0], ...
                         [p1, p2, p3], data)));
end
