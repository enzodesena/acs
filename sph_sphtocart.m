function x_data = sph_sphtocart(theta_phi_data)
% Converts from spherical coordinates to cartesian ones.
% Input should either be scalar or elements should be on a column vector
% and theta and phi should have the same length.
%    We are using a spherical coordinate system with the convention 
%    that theta is the angle formed with the z-axis, 
%    and phi is the angle formed on the projection on the x-y plane 
%    with the x-axis.
%    The x-axis corresponds to (r, theta, phi) = (r, pi/2, 0).

x_data = [sin(theta_phi_data(:, 1)).*cos(theta_phi_data(:, 2)),...
          sin(theta_phi_data(:, 1)).*sin(theta_phi_data(:, 2)),...
          cos(theta_phi_data(:, 1))];