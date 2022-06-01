%% optimize peak alignment for a multi-dynamic method

% for setup, use similar elements from Burdick's Beam Interpolation

% Iso_Name = {'Pb204','Pb205','Pb206','Pb207','Pb208'};
% F_ind = [0 0 0 0 0    2 3 4 5;...
%          0 0 0 0 1    3 4 5 0;...
%          0 0 0 1 2    4 5 0 0;...
%          0 0 1 2 3    5 0 0 0;...
%          0 1 2 3 4    0 0 0 0];

% KU dynamic Sr
Iso_Name = ["Sr84", "Rb85", "Sr86", "Sr87", "Sr88"];
Col_Name = ["L2", "L1", "H1", "H2"];
Iso_In_Axial = false;
F_massIndex = [1 2 3 4;
               0 3 4 5;
               0 4 5 0];
F_stationary = [false false true false];
massOffsetFromStationary = [-2 -1 0 1]; % can compute from F_ind if not lazy
input_AxMasses = [85.4; 86.4; 87.4];
Reff = 540;

% % Two-sequence dynamic Sr
% Iso_Name = {'Sr86', 'Sr87', 'Sr88'}; % keyed to class defined in mass.m
% Col_Name = {'Ax', 'H1'};             % collector names, # = cols of F_ind
% F_ind =    [1 2;...                  % indices correspond to position in
%             2 3];                    % Iso_Name
% Ax_Masses = [mass.Sr86; 
%              mass.Sr87];
% F_stationary = [true false];         % true if Faraday is stationary

% note: reference a mass as e.g., mass.Sr84 or mass.(Iso_Name{1})


%% set up problem

nIsotopes = length(Iso_Name); % # of isotopes measured
nCollectors = sum(any(F_massIndex)); % # of columns containing nonzero indices
nSequences = size(F_massIndex,1);

% temporary (?) assumption: all Iso_Name isotopes and collector (columns) 
% are used -- no unused isotopes or unused collectors.  also, one and only
% one stationary or declared-stationary Faraday. Positions solved relative
% to this stationary location

nPositions = nCollectors - 1; % see assumptions. # of unknowns


%% if no isotope in axial collector, sort out position of stationary one

if ~Iso_In_Axial

% first, find a sequence that includes the stationary collector
seqWithStationary = find(F_massIndex(:,F_stationary) ~= 0, 1);
massOnStationary = mass.(Iso_Name(F_massIndex(seqWithStationary,F_stationary)));

% how far away is the stationary collector from the axial position?
massInAxial = input_AxMasses(seqWithStationary);
massDifferenceStationaryAxial = massOnStationary - massInAxial;

% define the position of the stationary collector to reflext the mass
% difference
colPosStationary = massDifferenceStationaryAxial * ...
                   Reff / input_AxMasses(seqWithStationary);

else
    colPosStationary = 0; % Axial collector is at position 0
end % if ~Iso_In_Axial


%% set up optimization - new
% optimize for axial masses and collector positions together

nPeaks = sum(sum(F_massIndex > 0)); % number of equations/discrepancies to minimize
nVariables = nPositions + nSequences; % number of variables to solve for

% initialize positions of collectors
positionIndex = cumsum(~F_stationary).*~F_stationary;
pos0 = Reff/mean(input_AxMasses)*massOffsetFromStationary + colPosStationary;
pos0 = pos0(~F_stationary);

unk0 = [pos0'; input_AxMasses]; % initial guess
options = optimset('fminunc');
options.TolFun = 1e-10;
options.TolX = 1e-10;
options.MaxIter = 1e4;
%options.Display = 'iter';
%options.PlotFcns = @optimplotfval;

[soln, fval] = fminunc(@(unk) peakAlignmentPenalty(unk, nPositions, ... 
                nPeaks, F_massIndex, F_stationary, colPosStationary, ...
                Iso_Name, Reff), unk0, options);


%% visualize results




%% objective function

function sse = peakAlignmentPenalty(unk, nPositions, nPeaks, ...
                   F_massIndex, F_stationary, colPosStationary, ...
                   Iso_Name, Reff)

% parcel out unknowns in position and axial mass variables
unkPositions = unk(1:nPositions);
axMasses = unk(nPositions+1:end);
positions = 1:size(F_massIndex,2);
positions(F_stationary) = colPosStationary;
positions(~F_stationary) = unkPositions;

F_peakIndex = zeros(size(F_massIndex));
F_peakIndex(F_massIndex > 0) = 1:nPeaks; % index peaks in seq table

sse = 0;
for iPeak = 1:nPeaks

    [iSeq, jCol] = find(F_peakIndex == iPeak);
    iso_index = F_massIndex(iSeq, jCol);
    iMass = mass.(Iso_Name(iso_index));
    axMass = axMasses(iSeq);
    position = positions(jCol);
    discrep = (iMass - axMass) - position*axMass/Reff;

    sse = sse + discrep^2;


end % for iPeak

end % function sse
