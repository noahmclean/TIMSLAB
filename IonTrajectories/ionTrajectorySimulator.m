%% coordinates for extended geometry mass spec.
% units of mm, ion source at startPoint
% ions take off with angle theta from [1 0] (in degrees)

Reff = 540; % mm
boundaryRotation = 26.56; % degrees
boundaryAngleWithVertical = 45 - boundaryRotation;
R90 = [0 1; -1 0]; % rotation matrix

option = 3;
if option == 1 % three ions, three angles

    %startDirection = sqrt(2)/2*[1 1]';
    %radius = 540/2; % /2 for actual using extended geometry, 540 to start
    nTraces = 9;
    %radii = repmat([510 540 568.39]/2, 1, 3);
    radiusValues = [520 540 559]/2;
    radiusIndices = [1 2 3 1 2 3 1 2 3];
    radii = radiusValues(radiusIndices);
    %radii = repmat([530 540 549.5]/2, 1, 3);
    thetas = [43 45 47];
    startDirectionVectors = [cosd(thetas); sind(thetas)];
    startDirections = startDirectionVectors(:,[1 1 1 2 2 2 3 3 3]);
    cmap = [parula(5) 0.1*ones(5,1)];
    cmap = cmap(1:3,:);

elseif option == 2

    nTraces = 900; % divisible by 3
    radiusValues = [530 540 549.5]/2;
    radiusIndices = randi(3, [1, nTraces]);
    radii = radiusValues(radiusIndices);
    thetas = linspace(43, 47, nTraces/3); thetas = repmat(thetas, 1, 3);
    startDirections = [cosd(thetas); sind(thetas)];
    cmap = [parula(5) 0.1*ones(5,1)];
    cmap = cmap(1:3,:);

elseif option == 3

    nTracesPerIsotope = 50; % traces per isotope
    radiusValues = [510 520 530 540 550 560]/2;
    nIsotopes = length(radiusValues);
    nTraces = nTracesPerIsotope * nIsotopes;
    radiusIndices = randi(nIsotopes, [1, nTraces]);
    radii = radiusValues(radiusIndices);
    thetas = linspace(43, 47, nTracesPerIsotope); 
    thetas = repmat(thetas, 1, nIsotopes);
    startDirections = [cosd(thetas); sind(thetas)];

    %cmap = [parula(2*nIsotopes-1) 0.1*ones(2*nIsotopes-1,1)];
    %cmap = cmap(1:nIsotopes,:);
    cmap = [prism(nIsotopes) 0.1*ones(nIsotopes,1)];
end

%startPoint = sqrt( Reff^2 / 2)/2 * [1 1]';
startPoint = [0 0]';

% initialize drawing
figure('Position', [1 1 1200 800])
hA = axes();
set(hA, 'DataAspectRatio', [1 1 1], 'FontSize', 18)
xlim(hA, [0 1250]);
ylim(hA, [0 600]);
xlabel('length (mm)'), ylabel('length (mm)')
title('Extended Geometry, Effective Radius = 54 cm')
hold on


for iTrace = 1:nTraces

radius = radii(iTrace);

startDirection = startDirections(:,iTrace);

color = cmap(radiusIndices(iTrace),:);

%% determine left and right magnet boundaries


% left boundary
pointOnLBoundary = sqrt( Reff^2 / 2) * [1 1]';
lDirectionVector = [sind(boundaryAngleWithVertical) -cosd(boundaryAngleWithVertical)]';

% right boundary
pointOnRBoundary = pointOnLBoundary .* [2 1]';
rDirectionVector = [-sind(boundaryAngleWithVertical) -cosd(boundaryAngleWithVertical)]';


%% calculate intersection of ion free flight from source with L boundary

t = [startDirection -lDirectionVector]\(pointOnLBoundary - startPoint);
lIntersect = startPoint + startDirection*t(1);


%% calculate ciruclar path through magnetic field

% center of circle in circular path
% rotate path 90 deg to radius of circular path through magnet
directionToCircleCenter = R90*startDirection;
circleCenter = lIntersect + radius*directionToCircleCenter;

% intersection with left boundary (as theta)
thetaL = acos( (lIntersect(1) - circleCenter(1))/radius );

% intersection with right boundary (could be prettier, this from solver)
thetaR = -2*atan((radius*rDirectionVector(1) + (- pointOnRBoundary(1)^2*...
    rDirectionVector(2)^2 + 2*pointOnRBoundary(1)*pointOnRBoundary(2)*...
    rDirectionVector(1)*rDirectionVector(2) - 2*pointOnRBoundary(1)*...
    rDirectionVector(1)*rDirectionVector(2)*circleCenter(2) + 2*pointOnRBoundary(1)*...
    rDirectionVector(2)^2*circleCenter(1)-pointOnRBoundary(2)^2*rDirectionVector(1)^2+...
    2*pointOnRBoundary(2)*rDirectionVector(1)^2*circleCenter(2) - 2*...
    pointOnRBoundary(2)*rDirectionVector(1)*rDirectionVector(2)*circleCenter(1) + ...
    radius^2*rDirectionVector(1)^2 + radius^2*rDirectionVector(2)^2 - ...
    rDirectionVector(1)^2*circleCenter(2)^2 + 2*rDirectionVector(1)*...
    rDirectionVector(2)*circleCenter(1)*circleCenter(2) - rDirectionVector(2)^2*...
    circleCenter(1)^2)^(1/2))/(pointOnRBoundary(1)*rDirectionVector(2) - ...
    pointOnRBoundary(2)*rDirectionVector(1) + radius*rDirectionVector(2) - ...
    rDirectionVector(2)*circleCenter(1) + rDirectionVector(1)*circleCenter(2)));


%% free path from magnet to collector

rIntersect = circleCenter + radius * [cos(thetaR); sin(thetaR)];
colDirection = R90*(rIntersect-circleCenter)/norm(rIntersect-circleCenter);

% x-intercept
colPoint = rIntersect + colDirection * (-rIntersect(2)/colDirection(2));


%% make a drawing

% source
plot(hA, startPoint(1), startPoint(2), '.', 'MarkerSize', 20, 'MarkerEdgeColor', 'k')

% left boundary
magnetFaceHalfWidth = 50;
upperL = pointOnLBoundary - magnetFaceHalfWidth*lDirectionVector;
lowerL = pointOnLBoundary + magnetFaceHalfWidth*lDirectionVector;
xLeft = [upperL(1) lowerL(1)]; yLeft = [upperL(2) lowerL(2)];
line(hA, xLeft, yLeft, 'Color', 'r', 'LineWidth', 2)

% right boundary
upperR = pointOnRBoundary - magnetFaceHalfWidth*rDirectionVector;
lowerR = pointOnRBoundary + magnetFaceHalfWidth*rDirectionVector;
xRight = [upperR(1) lowerR(1)]; yRight = [upperR(2) lowerR(2)];
line(hA, xRight, yRight, 'Color', 'r', 'LineWidth', 2)

% upper and lower magnet boundaries, as circular arcs
magnetCenter = [sqrt( (Reff*3/2)^2/2 ); sqrt( (Reff/2)^2 / 2)];
centerToLowerLeft = magnetCenter - lowerL; centerToLowerRight = magnetCenter - lowerR;
centerToUpperLeft = magnetCenter - upperL; centerToUpperRight = magnetCenter - upperR;
magnetRadiusLower = norm(centerToLowerLeft);
magnetRadiusUpper = norm(centerToUpperLeft);
thetaLowerLeft  = acos(centerToLowerLeft(1) /magnetRadiusLower);
thetaLowerRight = acos(centerToLowerRight(1)/magnetRadiusLower);
thetaUpperLeft  = acos(centerToUpperLeft(1) /magnetRadiusUpper);
thetaUpperRight = acos(centerToUpperRight(1)/magnetRadiusUpper);
thetas = linspace(thetaLowerLeft, thetaLowerRight, 200);
plot(hA, magnetCenter(1) + magnetRadiusLower*cos(thetas), ...
         magnetCenter(2) + magnetRadiusLower*sin(thetas), 'Color', 'r', 'LineWidth', 2)
thetas = linspace(thetaUpperLeft, thetaUpperRight, 200);
plot(hA, magnetCenter(1) + magnetRadiusUpper*cos(thetas), ...
         magnetCenter(2) + magnetRadiusUpper*sin(thetas), 'Color', 'r', 'LineWidth', 2)

% ion free flight path from source
line(hA, [startPoint(1) lIntersect(1)], [startPoint(2) lIntersect(2)], 'Color', color, 'LineWidth', 1)

% ion path through magnet
thetas = linspace(thetaL, thetaR, 200);
plot(hA, circleCenter(1) + radius*cos(thetas), ...
     circleCenter(2) + radius*sin(thetas), 'Color', color, 'LineWidth', 1)

% ion path to collector
line(hA, [rIntersect(1) colPoint(1)], [rIntersect(2) colPoint(2)], 'Color', color, 'LineWidth', 1)


end % for iTrace
