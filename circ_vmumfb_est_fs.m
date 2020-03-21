function [opt_mu, opt_k1, opt_k2, opt_p1, opt_p2, opt_p3, max_ll] = ...
    circ_vmumfb_est_fs(data, mus, k1s, k2s, p1s, p2s)
%CIRC_VMUM_EST_FS full search estimation of vMUM model parameters
% 
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena


assert(iscolumn(data));

max_ll = - Inf;
opt_mu = nan;
opt_k1 = nan;
opt_k2 = nan;
opt_p1 = nan;
opt_p2 = nan;
opt_p3 = nan;

for mus_i = 1:length(mus)
    for k1_i = 1:length(k1s)
        for k2_i = 1:length(k2s)
            for p1s_i = 1:length(p1s)
                for p2s_i = 1:length(p2s)
                    for p3s_i = 1:length(p2s)
                        p3 = 1 - p1s(p1s_i) - p2s(p2s_i);

                        if p3 < 0 || p3 > 1
                            continue
                        end
                        %% Calculate
                        ll = circ_vmumfb_ll(mus(mus_i), k1s(k1_i), k2s(k2_i), p1s(p1s_i), p2s(p2s_i), p3, data, false, false);
                        
                        if ll > max_ll
                            max_ll = ll;
                            opt_mu = mus(mus_i);
                            opt_k1 = k1s(k1_i);
                            opt_k2 = k2s(k2_i);
                            opt_p1 = p1s(p1s_i);
                            opt_p2 = p2s(p2s_i);
                            opt_p3 = p3;
                        end
                    end
                end
            end    
        end
    end
end

