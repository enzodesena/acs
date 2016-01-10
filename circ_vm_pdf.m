function pdf = circ_vm_pdf(mu, k, theta)
%CIRC_VM_PDF von Mises distribution PDF
%   pdf = CIRC_VM_PDF(mu, k, theta) returns the values of the PDF
%   of the von Mises distribution with parameters mu and k
%   at angle(s) theta
%   
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

pdf = 1/(2*pi*besseli(0,k))*exp(k*cos(theta-mu));
