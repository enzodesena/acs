function sph_hist(theta_phi_data, nbins)

if nargin == 1
    nbins = 30;
end

thetas = linspace(0, pi, nbins);
phis = linspace(0, 2*pi, nbins);

[~, D] = size(theta_phi_data);
assert(D == 2);

N = hist3(theta_phi_data, {thetas, phis});
size(N)

[f, t] = meshgrid(phis, thetas);
x = sin(t)*cos(f);
y = sin(t)*sin(f);
z = cos(t);
subplot(1,2,1)
surf(x,y,z, N, 'edgecolor','none');

subplot(1,2,2)
surf(f,t,N, 'edgecolor','none');