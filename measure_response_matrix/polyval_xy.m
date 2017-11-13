% [z] = polyval_xy(x, y, coeffs, m)
% 
% Calculates "z"("x","y")for a bivariate polynomial with coefficients
% defined by "coeffs" and orders "m", as determined during the DM 
% calibration porocess.  

function [z] = polyval_xy(x, y, coeffs, m)

% Find size of input arrays; size x and y must be the same
[rows,cols] = size(x);
z = zeros(rows,cols);

% Loop through orders to compute poly val
for j = 1:length(m)
    z = z + coeffs(j)*x.^m(j,1).*y.^m(j,2);
end

