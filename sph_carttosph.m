function [theta_phi_data] = sph_carttosph(x)
% Converts from spherical coordinates to cartesian ones.
% Input should either be scalar or elements should be on a column vector.
%    We are using a spherical coordinate system with the convention 
%    that theta is the angle formed with the z-axis, 
%    and phi is the angle formed on the projection on the x-y plane 
%    with the x-axis.
%    The x-axis corresponds to (r, theta, phi) = (r, pi/2, 0).

theta_phi_data = [acos(x(:, 3)), atan(x(:, 2)./x(:, 1))];