function [mu_hat, k_hat, p1_hat, p2_hat] = circ_vmum_est_mm(data)

phi = mod(data*2, 2*pi);

%% Estimate mu
mu_hat = circ_mean(phi)/2;

%% Estimate k
ac2w = mean(cos(2*(phi-2*mu_hat)));
ac1w = mean(cos(1*(phi-2*mu_hat)));


%options = optimoptions('fsolve','TolFun',1E-9);
k_hat = fsolve(@(k)besseli(4,k)./besseli(2,k)*ac1w-ac2w, 1);


%% Assumes prior knowledge that k<100
%k_hat = fsolve(@(k)besseli(4,k)./besseli(2,k)-min(abs(ac2w/ac1w), 0.94), 1, options);

%% Works worse because rho2w>ac2w
% bc2w = mean(sin(2*(phi-2*mu_hat))); rho2w = sqrt(ac2w.^2+bc2w.^2);
% k_hat = fsolve(@(k)ac1w.*besseli(4,k)./besseli(2,k)-rho2w, 1, options);

%k_hat = fsolve(@(k)ac1w*besseli(4,k)./besseli(2,k)-ac2w, 1, options); 

%% Version that aims to remove bias
%k_hat = fsolve(@(k)ac1w*besseli(4,k)./besseli(2,k)-ac2w, 1)-4.54E-5*exp(15.06*abs(ac2w/ac1w)); 

%% Doesn't work because besseli(4,k) and besseli(2,k) are very large numbers
%k_hat = fsolve(@(k)abs(ac1w)*besseli(4,k)-besseli(2,k)*abs(ac2w), 1); 

p_ol = min(abs(ac1w)*besseli(0,k_hat)/besseli(2,k_hat), 1);

%% Sequential correction
k_hat = fsolve(@(k)besseli(2,k)./besseli(0,k)*p_ol-(abs(ac1w)), 1);


%% Works worse than other methods
% [mean(cos(1*(phi-2*mu_hat))),mean(cos(2*(phi-2*mu_hat))),mean(cos(3*(phi-2*mu_hat))),mean(cos(4*(phi-2*mu_hat)))]
% [mean(sin(1*(phi-2*mu_hat))),mean(sin(2*(phi-2*mu_hat))),mean(sin(3*(phi-2*mu_hat))),mean(sin(4*(phi-2*mu_hat)))]
% params = fmincon(@(params)(...
%     (besseli(2,params(1))./besseli(0,params(1))*params(2)-mean(cos(1*(phi-2*mu_hat)))).^2+...
%     (besseli(4,params(1))./besseli(0,params(1))*params(2)-mean(cos(2*(phi-2*mu_hat)))).^2+...
%     (besseli(6,params(1))./besseli(0,params(1))*params(2)-mean(cos(3*(phi-2*mu_hat)))).^2+...
%     (besseli(8,params(1))./besseli(0,params(1))*params(2)-mean(cos(4*(phi-2*mu_hat)))).^2), [1, 0.5],...
%     [],[],[],[],[0, 0],[inf,1]);
% 
% k_hat = params(1);
% p_ol = params(2);

%% Estimate p1 and p2
a1 = mean(cos(data-mu_hat));
p1_hat = min(abs(p_ol + a1*besseli(0,k_hat)/besseli(1,k_hat))/2, p_ol);
p2_hat = p_ol - p1_hat;

%% Asserts
assert(p1_hat>=0 && p1_hat<=1);
assert(p2_hat>=0 && p2_hat<=1);
assert((1-p1_hat-p2_hat)>=0 && (1-p1_hat-p2_hat)<=1);


