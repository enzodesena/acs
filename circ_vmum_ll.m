function ll = circ_vmum_ll(mu, k, p1, p2, p3, data, relaxed)
%CIRC_VMUM_LL returs log-likelihood of the vMUM model
%
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

%% Asserts
assert(isscalar(mu));
assert(isscalar(k));
assert(isscalar(p1));
assert(isscalar(p2));
assert(isscalar(p3));
assert(isvector(data));

if nargin == 7 && relaxed 
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
    circ_vmum_asserts(mu, k, p1, p2, p3, 1e-3);
end


%% Calculate
ll = sum(log(circ_vmum_pdf(mu, k, p1, p2, p3, data)));
assert(not(isnan(ll)));

end
