function [mu_hat, k_hat, p1_hat, p2_hat, p3_hat, ll] = circ_vmum_est_ml(data, x0)
options = optimoptions('fmincon', 'Display', 'off');

%% Run
params = fmincon(@(params) ...
           -getll(params, data), ...
           x0, [], [], [0, 0, 1, 1, 1], 1, ...
           [-pi, 0, 0, 0, 0], [pi, inf, 1, 1, 1], [], options);
       
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