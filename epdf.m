function [pdf_thetas, pdf_data] = epdf(data, values)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[N, pdf_thetas] = hist(data, values);
pdf_data =  N./trapz(pdf_thetas, N);

if nargout == 0
    plot(pdf_thetas, pdf_data)
end

end

