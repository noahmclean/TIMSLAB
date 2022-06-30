function B = bbase(x, xl, xr, nseg, bdeg)
%BBASE Assemble B-spline basis matrix with evenly spaced knots
%   adapted from bbase.R in CRAN package JOPS
%
%   x is the vector of argument values at which the B-spline basis vectors
%   are to be evaluated
%   xl is the lower limit of the domain of x
%   xr is the upper limit on the domain of x
%   nseg is the number of equally spaced segments
%   bdeg is the degree of the spline (3 is cubic)
%
%   B is the B-spline basis matrix, such that u = B*x
%   B has length(x) rows and nseg+bdeg columns
%   R authors: Paul Eilers and Brian Marx

% Check domain, adjust if necessary
if xl > min(x)
    xl = min(x);
    disp('Adjusted left boudary to min(x)')
end
if xr < max(x)
    xr = max(x);
    disp('Adjusted right boundary to max(x)')
end

% make x a column vector
x = x(:);

% Compute the B-splines
dx = (xr-xl)/nseg;
knots = linspace(xl-bdeg*dx, xr+bdeg*dx, nseg+2*bdeg+1); % xl-bdeg*dx:dx:x4+bdeg*dx;

nx = length(x);
nt = length(knots);
X = kron(x, ones(1,nt));
T = kron(knots, ones(nx,1));

P = (X - T).^bdeg .* (X >= T); % added = to >= to match de Boor definition?
D = diff(eye(nt), bdeg+1)/(gamma(bdeg+1)*dx^bdeg);
B = (-1)^(bdeg+1) * P * D';

% Make B-splines exactly zero beyond their end knots
nb = size(B,2);
X = kron(x, ones(1,nb));
sk = knots((1:nb) + bdeg + 1);
SK = kron(sk, ones(nx,1));
Mask = X < SK;
B = B .* Mask;

end % function bbase

