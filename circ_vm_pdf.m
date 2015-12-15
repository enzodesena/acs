function pdf = circ_vm_pdf(mu, k, theta)

pdf = 1/(2*pi*besseli(0,k))*exp(k*cos(theta-mu));
