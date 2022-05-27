%% optimize peak alignment for a multi-dynamic method

% for setup, use similar elements from Burdick's Beam Interpolation

% Iso_Name = {'Pb204','Pb205','Pb206','Pb207','Pb208'};
% F_ind = [0 0 0 0 0    2 3 4 5;...
%          0 0 0 0 1    3 4 5 0;...
%          0 0 0 1 2    4 5 0 0;...
%          0 0 1 2 3    5 0 0 0;...
%          0 1 2 3 4    0 0 0 0];
%
% KU dynamic Sr
Iso_Name = ["Sr84", "Rb85", "Sr86", "Sr87", "Sr88"];
Col_Name = ["L2", "L1", "H1", "H2"];
Iso_In_Axial = false;
F_ind = [1 2 3 4;
         0 3 4 5;
         0 4 5 0];
F_stationary = [false false true false];
Ax_Masses = [85.4; 86.4; 87.4];
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
nCollectors = sum(any(F_ind)); % # of columns containing nonzero indices
nSequences = size(F_ind,1);

% temporary (?) assumption: all Iso_Name isotopes and collector (columns) 
% are used -- no unused isotopes or unused collectors.  also, one and only
% one stationary or declared-stationary Faraday.


nPositions = nCollectors - 1; % see assumptions. # of unknowns


%% determine number of equations

Unk_ind = zeros(size(F_ind));
nEquations = 0;
for iSeq = 1:nSequences
    for jCol = 1:nCollectors

    if F_ind(iSeq,jCol) > 0 && ~F_stationary(jCol)
        nEquations = nEquations + 1;
        Unk_ind(iSeq,jCol) = nEquations;

    end % if

    end
end

%% if no isotope in axial collector, position stationary collector 

% first, find a sequence that includes the stationary collector
seqWithStationary = find(F_ind(:,F_stationary) ~= 0, 1);
massOnStationary = mass.(Iso_Name(F_ind(seqWithStationary,F_stationary)));

% how far away is the stationary collector from the axial position?
massInAxial = Ax_Masses(seqWithStationary);
massDifferenceStationaryAxial = massOnStationary - massInAxial;

% define the position of the stationary collector to reflext the mass
% difference
colPosStationary = massDifferenceStationaryAxial * ...
                                          Reff/Ax_Masses(seqWithStationary);


% adjust the axial masses for other sequences to reflect cent 


%% set up optimization

% b is the (amu) distance between axial and measured masses
b = zeros(nEquations,1);
% A relates position of the cup (in mm) to peak center for a given mass
A = zeros(nEquations, nCollectors);
% colPos = A*b, where colPos is distances from axial for each cup
for iEqn = 1:nEquations

    [iSeq, jCol] = find(Unk_ind == iEqn);
    iso_index = F_ind(iSeq, jCol);
    iMass = mass.(Iso_Name(iso_index));
    AxMass = Ax_Masses(iSeq);
    b(iEqn) = iMass - AxMass;
    A(iEqn, jCol) = AxMass/Reff;

end % for iEqn

% clear column of A assigned to stationary detector
A = A(:, ~F_stationary);
colPos = A\b; % collector positions, in mm from axial position


function 