function [opt_mu, opt_ks, opt_p1, opt_p2, opt_p3, max_ll] = ...
    circ_vmum_est_fs(data, mus, ks, p1s, p2s)
%CIRC_VMUM_EST_FS full search estimation of vMUM model parameters
% 
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena


assert(iscolumn(data));

max_ll = - Inf;
opt_mu = nan;
opt_ks = nan;
opt_p1 = nan;
opt_p2 = nan;
opt_p3 = nan;

for mus_i = 1:length(mus)
    for ks_i = 1:length(ks)
        for p1s_i = 1:length(p1s)
            for p2s_i = 1:length(p2s)
                for p3s_i = 1:length(p2s)
                    p3 = 1 - p1s(p1s_i) - p2s(p2s_i);
                    
                    if p3 < 0 || p3 > 1
                        continue
                    end
                    %ll(mus_i, ks_i, p1s_i, p2s_i, p3s_i)
                    ll =  sum(log(circ_vmm_pdf([mus(mus_i), mus(mus_i)+pi, 0], ...
                                              [ks(ks_i), ks(ks_i), 0], ...
                                              [p1s(p1s_i), p2s(p2s_i), p3], data)));
                    if ll > max_ll
                        max_ll = ll;
                        opt_mu = mus(mus_i);
                        opt_ks = ks(ks_i);
                        opt_p1 = p1s(p1s_i);
                        opt_p2 = p2s(p2s_i);
                        opt_p3 = p3;
                    end
                end
            end
        end
    end
end

