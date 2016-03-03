function [pdf_thetas, pdf_data] = circ_hist(data, resolution, linespec)
%CIRC_HIST circular histogram
%   [pdf_thetas, pdf_data] = CIRC_HIST(data, resolution) returns
%   the x axis values pdf_thetas and y axis values for a histogram
%   of given data and resolution
% 
%   Audio Circular Statistics (ACS) library
%   Copyright 2016 Enzo De Sena

if nargin == 1
    resolution = pi/100;
end

if nargin < 2
    linespec = '';
end

assert(iscolumn(data));
assert(isscalar(resolution));

data = wrappi(data);
[N, pdf_thetas] = hist(data, -pi:resolution:pi);

pdf_data =  N./trapz(pdf_thetas, N);
pdf_data(1) = pdf_data(1) * 2; % The first and last bin have half size
pdf_data(end) = pdf_data(end) * 2;

if nargout < 2
    plot(pdf_thetas, pdf_data, linespec);
end

end
