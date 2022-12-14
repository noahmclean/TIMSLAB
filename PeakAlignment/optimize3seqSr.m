%% optimize three sequence Sr method as a test/sandbox

%     L3  L2  Ax  H1  H2
% S1          86  87  88
% S2      86  87  88
% S3  86  87  88

% parameterized as 
% m1, m2, m3 axial masses during S1, S2, S3
% d1, d2, d3, d4 distances between L2-L2, L2-Ax, Ax-H1, H1-H2

% constants to use: a,b,c are inverse of dispersion for S1, S2, S3
a = 85.9121/540;
b = 86.9031/540;
c = 87.9085/540;

y = repmat([mass.Sr86; mass.Sr87; mass.Sr88], 3, 1);
A = [1 0 0  0  0  0  0;
     1 0 0  0  0  a  0;
     1 0 0  0  0  a  a;
     0 1 0  0 -b  0  0;
     0 1 0  0  0  0  0;
     0 1 0  0  0  b  0;
     0 0 1 -c -c  0  0;
     0 0 1  0 -c  0  0;
     0 0 1  0  0  0  0];
x = A\y;
r = y - A*x;

% results in terms of peak top width
WC = 1.0; % mm, width of collector slit
WB = 0.35; % mm, width of beam with double focusing
Reff = 540; % mm, effective radius of magnet
WT = (WC - WB)*[mass.Sr86 mass.Sr87 mass.Sr88]/540;