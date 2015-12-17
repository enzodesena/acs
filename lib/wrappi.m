
% wraps to -pi pi
function wrapped_angle = wrapToMinusPiPi(angle)
wrapped_angle = mod(angle+pi,2*pi)-pi;
end