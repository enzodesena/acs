% wraps to -pi pi
function wrapped_angle = wrappi(angle)
wrapped_angle = rem(angle, 2*pi);
above_pi = abs(wrapped_angle) > pi;
wrapped_angle(above_pi) = wrapped_angle(above_pi) - 2*pi*sign(wrapped_angle(above_pi));
end