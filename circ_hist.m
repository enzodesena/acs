function [pdf_thetas, pdf_data] = circ_hist(data, resolution, linespec)

if nargin == 2
    linespec = '';
end

%data = wrapToMinusPiPi(data(:));
[N, pdf_thetas] = hist(data, -pi:resolution:pi);

pdf_data =  N./trapz(pdf_thetas, N);
pdf_data(1) = pdf_data(1) * 2; % The first and last bin have half size
pdf_data(end) = pdf_data(end) * 2;

if nargout < 2
    plot(pdf_thetas, pdf_data, linespec);
end

end

% wraps to -pi pi
function wrapped_angle = wrapToMinusPiPi(angle)
wrapped_angle = mod(angle+pi,2*pi)-pi;
end