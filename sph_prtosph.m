function theta_phi_data = sph_prtosph(pitch_roll_data)

[N, D] = size(pitch_roll_data);
assert(D == 2);

pitch = pitch_roll_data(:, 1);
roll = pitch_roll_data(:, 2);

theta_phi_data = nan(N, 2);

for n=1:N
    theta = mod(pitch(n), 2*pi);
    phi = mod(roll(n), 2*pi);
    if theta > pi
        phi = phi+pi;
        theta = 2*pi-theta;
    end
    phi = mod(phi, 2*pi);
    
    assert(theta>0 & theta<=pi);
    assert(phi>0 & phi<=2*pi);
    
    theta_phi_data(n, :) = [theta, phi];
end
