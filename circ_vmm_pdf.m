function pdf = circ_vmm_pdf(mu, k, p, theta)

pdf = zeros(size(theta));
for i=1:length(mu)
    pdf = pdf + p(i)/(2*pi*besseli(0,k(i)))*exp(k(i)*cos(theta-mu(i)));
end
