function [front_data, back_data, err_data, front_idxs, back_idxs, err_idxs] = ...
    circ_vmum_comp(data)
%CIRC_VMUM_COMP returns an estimate of the three groups of data in the 
%   vMUM distibution, that is, frontal images, back images and 
%   uniformly distributed errors.
% 
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

assert(iscolumn(data));

[mu_hat, k_hat, p1_hat, p2_hat, p3_hat] = circ_vmum_est_ml(data);

l_front = p1_hat*circ_vm_pdf(mu_hat, k_hat, data);
l_back = p2_hat*circ_vm_pdf(mu_hat+pi, k_hat, data);
l_err = p3_hat/(2*pi);

front_idxs = l_front >= l_back & l_front >= l_err;
front_data = data(front_idxs);
back_idxs = l_back >= l_front & l_back >= l_err;
back_data = data(back_idxs);
err_idxs = l_err > l_front & l_err > l_back;
err_data = data(err_idxs);

end