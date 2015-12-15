function pdf = circ_vmm3_pdf(mu_hat, k_hat, p, theta)

pdf = 1/(2*pi*besseli(0,k_hat))*(p*exp(k_hat*cos(theta-mu_hat))+...
    (1-p)*exp(-k_hat*cos(theta-mu_hat)));
