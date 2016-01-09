function [mu_hat, k_hat, p1_hat, p2_hat, p3_hat, ll] = ...
    circ_vmum_est_ml(data, mu_0, k_0, p1_0, p2_0, p3_0, verb)
if nargin == 6 || verb == false
    options = optimoptions('fmincon', 'Display', 'off');
else
    options = optimoptions('fmincon');
end

%% Asserts
assert(p1_0>=0.0 & p1_0<=1.0);
assert(p2_0>=0.0 & p2_0<=1.0);
assert(p3_0>=0.0 & p3_0<=1.0);
assert(abs((p1_0+p2_0+p3_0)-1.0)<1E-10);
assert(isscalar(mu_0) & isscalar(k_0));
assert(isscalar(p1_0) & isscalar(p2_0) & isscalar(p3_0));
assert(isvector(data));

%% Adjust initial point to problem domain
if k_0 < 0
    k_0 = -k_0;
    mu_0 = mu_0 + pi;
    temp = p1_0;
    p1_0 = p2_0;
    p2_0 = temp;
end
mu_0 = mod(mu_0, 2*pi);

%% Run
params = fmincon(@(params) ...
           -getll(params, data), ...
           [mu_0, k_0, p1_0, p2_0, p3_0], [], [], [0, 0, 1, 1, 1], 1, ...
           [0, 0, 0, 0, 0], [2*pi, inf, 1, 1, 1], [], options);
       
%% Assign output values
ll = getll(params, data);

mu_hat = params(1);
k_hat  = params(2);
p1_hat = params(3);
p2_hat = params(4);
p3_hat = params(5);

end

function l = getll(params, data)
l = sum(log(circ_vmm_pdf([params(1), params(1)+pi, 0], ...
                         [params(2), params(2), 0], ...
                         [params(3), params(4), params(5)], data)));
end