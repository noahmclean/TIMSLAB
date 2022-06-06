%% coordinates for extended geometry mass spec.
% units of mm, ion source at (0,0), 
% ions take off to NE on line y = x

Reff = 540; % mm
boundaryRotation = 26.56; % degrees
boundaryAngleWithVertical = 45 - boundaryRotation;

startDirection = sqrt(2)/2*[1 1]';
radius = 540/2; % /2 for actual using extended geometry, 540 to start

R90 = [0 1; -1 0]; Rminus90 = [0 -1; 1 0]; % rotation matrices


%% determine left and right magnet boundaries


% left boundary
pointOnLBoundary = sqrt( Reff^2 / 2) * [1 1]';
lDirectionVector = [sind(boundaryAngleWithVertical) -cosd(boundaryAngleWithVertical)]';

% right boundary
pointOnRBoundary = pointOnLBoundary .* [2 1]';
rDirectionVector = [-sind(boundaryAngleWithVertical) -cosd(boundaryAngleWithVertical)]';


%% calculate intersection of ion free flight from source with L boundary

t = [startDirection -lDirectionVector]\pointOnLBoundary;
lIntersect = startDirection*t(1);


%% calculate ciruclar path through magnetic field

% center of circle in circular path
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

rIntersect = 


%% make a drawing

figure('Position', [1 1 1200 900])
hA = axes();
set(hA, 'DataAspectRatio', [1 1 1])
xlim(hA, [0 1150]);
ylim(hA, [0 800]);
hold on

% source
plot(0, 0, '.', 'MarkerSize', 20, 'MarkerEdgeColor', 'k')
%focus
plot(sqrt((3*Reff)^2 /2), 0, '.', 'MarkerSize', 20, 'MarkerEdgeColor', 'r')

% left boundary
magnetFaceHalfWidth = 50;
upperL = pointOnLBoundary + magnetFaceHalfWidth*lDirectionVector;
lowerL = pointOnLBoundary - magnetFaceHalfWidth*lDirectionVector;
xLeft = [upperL(1) lowerL(1)]; yLeft = [upperL(2) lowerL(2)];
line(xLeft, yLeft, 'Color', 'r', 'LineWidth', 2)

% right boundary
upperR = pointOnRBoundary + magnetFaceHalfWidth*rDirectionVector;
lowerR = pointOnRBoundary - magnetFaceHalfWidth*rDirectionVector;
xRight = [upperR(1) lowerR(1)]; yRight = [upperR(2) lowerR(2)];
line(xRight, yRight, 'Color', 'r', 'LineWidth', 2)

% ion free flight path from source
line([0 lIntersect(1)], [0 lIntersect(2)], 'Color', 'k', 'LineWidth', 1)

% ion path through magnet
circleStep = -pi/100;
thetas = (thetaL:circleStep:thetaR)';
plot(circleCenter(1) + radius*cos(thetas), ...
     circleCenter(2) + radius*sin(thetas), 'Color', 'k', 'LineWidth', 1)