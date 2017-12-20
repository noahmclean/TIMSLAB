%% Interpret multiple Pb standards run on the Phoenix

% Pb standard runs should be saved as a .xlsx folder and all saved in one folder.
% That folder should be dataFileFolder.
% The name of each .xlsx data file should start with NBS981 or NBS982.

dataFileFolder = '/Users/noahmc/Documents/KU/IGL/PhoenixData/PhoenixPbStandardData_SEM02/';

MSmethod.measMasses = {'204', '205', '206', '207', '208'}; 
MSmethod.outRatios = {'204/206', '205/206', '207/206', '208/206'}; 
MSmethod.cyclesPerBlock = 12; 
MSmethod.BItimes982 = [10 2 1 2 3 2 5 2 3 2];
MSmethod.BItimes981 = [10 2 1 2 5 2 5 2 2 2];
MSmethod.BItimes = MSmethod.BItimes982;
MSmethod.BImethod = 'Dodson'; %options: 'Dodson', 'Quadrift'


% Parse folder of Pb standard data files
runs = parsePbStandardsPhoenix(dataFileFolder);
n.runs = size(runs,1);

% Do beam interpolation
for irun = 1:n.runs
    
    switch runs(irun).standard
        case 'NBS981'
            MSmethod.BItimes = MSmethod.BItimes981;
        case 'NBS982'
            MSmethod.BItimes = MSmethod.BItimes982;
    end % switch
    
    switch MSmethod.BImethod
        case 'Dodson'
            runs(irun).BIdata = DodsonBI(runs(irun).dataDT0, MSmethod);
        case 'Quadrift'
            runs(irun).BIdata = QuadDriftCorr_v1(runs(irun).dataDT0, MSmethod);
    end % switch
    
end % for

% Plot data and do brushing
runs = brushPbStandardData_v1(runs, MSmethod);

% brushPbStandardData will handoff to reduceBrushedPbData after UI button press