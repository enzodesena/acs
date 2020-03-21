function ll = circ_vmumfb_ll(mu, k1, k2, p1, p2, p3, data, check_values, relaxed)
%CIRC_VMUMFB_LL returs log-likelihood of the vMUMFB model
%
%   Audio Circular Statistics (ACS) library
%   Copyright 2020 Enzo De Sena

%% Asserts
if nargin <= 7 || check_values
    assert(isscalar(mu));
    assert(isscalar(k1));
    assert(isscalar(k2));
    assert(isscalar(p1));
    assert(isscalar(p2));
    assert(isscalar(p3));
    assert(isvector(data));

    if nargin == 9 && relaxed 
        if abs(p1+p2+p3 - 1.0)>1e-3
            warning(strcat('Probabilities where normalised because they ',...
                'did not sum to one: ', ...
                num2str(p1), ' ', num2str(p2), ' ', num2str(p3), ' sum: ', ...
                num2str(p1+p2+p3)));
            temp = p1+p2+p3;
            p1 = p1/temp;
            p2 = p2/temp;
            p3 = p3/temp;
        end
    else
        circ_vmumfb_asserts(mu, k1, k2, p1, p2, p3, 1e-3);
    end
end
 
%% Calculate
ll = sum(log(circ_vmm_pdf([mu, pi-mu, 0], [k1, k2, 0], [p1, p2, p3], data)));
assert(not(isnan(ll)));

end
