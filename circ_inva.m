function y = circ_inva(R)
%CIRC_INVA Approximation of inverse Bessel
%   y = CIRC_INVA(R) Returns an approximation of the solution to the 
%   problem besseli(1, x)/besseli(0,x) - R = 0
%   using approximation given in (5.3.11)
%   in Mardia, "Directional Statistics"
% 
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

assert(all(R >= 0 & R <= 1));
y = (1.28 - 0.53.*R.^2).*tan(pi.*R./2);

% Banerjee et al.
% "Clustering on the Unit Hypersphere using von Mises-Fisher Distributions"
% eq (4.4)
%y = (R.*2-R.^3)./(1-R.^2);

% Approx in Mardia, "Directional Statistics"
% eq (5.3.7), (5.3.8) and (5.3.10)
% % y = (R < 0.53) .* (2.*R + R.^3 + 5/6.*R.^5) + ...
% %     (R >= 0.53) .* (R< 0.85) .* (1./(2.*(1-R)-(1-R).^3-(1-R).^3)) + ...
% %     (R >= 0.85) .* (-0.4+1.39.*R+0.43./(1-R));
end