classdef splineBasisModel
    %SPLINEBASISMODEL B-spline basis
    %   using bbase() from R Package from Joy of P-Splines
    %   R authors: Paul Eilers and Brian Marx
    
    properties (Access = protected)
        bdeg % degree of basis function (1 is linear, 2 is quadratic, 3 is cubic)
        nseg % number of equally spaced segments on which to fit x variable
    end

    properties (SetAccess = protected)
        beamMassInterp  % vector of x values at which to evaluate spline
        B               % B-spline basis matrix, with length(x) rows, nseg+bdeg columns
    end
    
    methods
        function splineBasis = splineBasisModel(beamMassInterp, nseg, bdeg)
            %SPLINEBASISMODEL Construct an instance of this class
            %   Detailed explanation goes here
            
            arguments
                beamMassInterp {mustBeVector}
                nseg (1,1) double = 1000
                bdeg (1,1) double = 3
            end

            splineBasis.beamMassInterp = beamMassInterp;
            splineBasis.nseg = nseg;
            splineBasis.bdeg = bdeg;
            splineBasis.B = bbase(beamMassInterp, nseg, bdeg);

        end % constructor function
        
    end % methods
end % classdef



function B = bbase(x, nseg, bdeg)
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

%   Note: elminated xl and xr arguments for peak center project
%   Added bbase() to splineBasisModel

% make x a column vector
x = x(:);
xl = x(1); 
xr = x(end);

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



