%% optimize three sequence Sr method as a test/sandbox

%     L3  L2  Ax  H1  H2
% S1          86  87  88
% S2      86  87  88
% S3  86  87  88

% parameterized as 
% m1, m2, m3 axial masses during S1, S2, S3
% d1, d2, d3, d4 distances between L2-L2, L2-Ax, Ax-H1, H1-H2

% constants to use: a,b,c are inverse of dispersion for S1, S2, S3
% note adjusted by ~0.0001-2 for offsets m1, m2, m3
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

isotopeMasses = [mass.Sr86 mass.Sr87 mass.Sr88];
collectorWidths = [0.87 0.87 0.87 0.87 0.87];

%% reproduce a mass scan using this information

% start simple with a beam shape from beamShape in TIMSLAB
% need: beamDistInterp: x-axis for beam profile, in mm distances
% need: beamShape: fit by beamShape routine

nMassesToScan = 700;
massScan = linspace(85.7, 86.3, nMassesToScan);

collectorPositions = [x(4)+x(5) x(4) 0 x(6) x(7)];
collectorOpenings = [collectorPositions - 0.5*collectorWidths;
                     collectorPositions + 0.5*collectorWidths];


for iMass = 1:nMassesToScan

    axialMass = massScan(iMass);
    dispersion = Reff/axialMass;
    isotopePositions = (isotopeMasses - axialMass)*dispersion;



end % for iMass
