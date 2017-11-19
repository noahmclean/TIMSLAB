%% Interpret multiple Pb standards run on the Phoenix

% Pb standard runs should be saved as a .xlsx folder and all saved in one folder.
% That folder should be inside the folder containing this script.
% The name of each .xlsx data file should start with NBS981 or NBS982.

dataFileFolder = 'PhoenixPbStandardData/'; % note, make sure to end string with /
BImethod = 'Dodson'; %options: 'Dodson', 'Quadrift'

MSmethod.measMasses = {'204', '205', '206', '207', '208'}; 
MSmethod.outRatios = {'204/206', '206/205', '207/206', '208/206'}; 
MSmethod.cyclesPerBlock = 12; 
MSmethod.BItimes = [10 2 1 2 3 2 5 2 3 2];

% Parse folder of Pb standard data files
runs = parsePbStandardsPhoenix(dataFileFolder);
n.runs = size(runs,1);

% Do beam interpolation
for ii = 1:n.runs
    switch BImethod
        case 'Dodson'
            runs(ii).BIdata = DodsonBI_v1(runs(ii).dataDT0, MSmethod);
        case 'Quadrift'
            runs(ii).BIdata = DodsonBI_v1(runs(ii).dataDT0, MSmethod);
    end % switch
end % for

% Plot data and do brushing
runs = brushPbStandardData_v1(runs, MSmethod);

