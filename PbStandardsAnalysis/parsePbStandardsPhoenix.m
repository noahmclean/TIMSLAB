function runs = parsePbStandardsPhoenix(folderstring)

%% parse a folder of 981 runs from the x62

% create the structure 'runs', with all relevant run data and parameters
% Note: parsei structure contains field 'isosHere' with isotopes present in each output function: USEFUL!

%% create the data structure

% write folderstring as 'folderName/'
fileStruct = dir([folderstring '*.xlsx']); %make a structure with properties of all .xlsx files in folder
fileN = numel(fileStruct); %this is the number of .xlsx files

runs = struct('name', {}, 'dataRaw', {}, 'dt', {}, 'time', {}, 'bi', {}, 'dataDT0', {}, 'secs', {}, 'temp', {});

for i = 1:fileN
    
    disp(['reading ', fileStruct(i).name, ', file no. ' , num2str(i)])
    %import cycle data
    [temp.numdata, temp.strdata] = xlsread( [folderstring fileStruct(i).name], 'CYCLE'); %cycle data
    
    %figure out function names (as cell array)
    temp.numfun = xlsread( [folderstring fileStruct(i).name], 'CTRL', 'B13' ); %number of functions
    temp.cycleColStart = 3; % first column of cycle data in CYCLE worksheet
    temp.cycleColEnd   = temp.cycleColStart + temp.numfun - 1;
    %temp.funstr = ['C9:C' num2str(9+temp.numfun-1)]; %cells of function names in IsoWorks sheet
    parsei.nameCell = temp.strdata(2, temp.cycleColStart : temp.cycleColEnd);
    %[~, parsei.nameCell] = xlsread( ['xlsxs/' fileStruct(i).name], 'IsoWorks', temp.funstr ); %function names
    
    % Read in following datafile properties:
    %name:
    [~, runs(i).name] = xlsread( [folderstring fileStruct(i).name], 'CTRL', 'D63' );
    runs(i).name = cell2mat(runs(i).name);
    %dead time:
    runs(i).dt   = xlsread( [folderstring fileStruct(i).name], 'TUNING', 'J10' );
    %timestamp, converted to MATLAB's datenum:
    [~, runs(i).time] = xlsread( [folderstring fileStruct(i).name], 'CTRL', 'D22' );
    runs(i).time = cell2mat(runs(i).time);
    temp.spaces = strfind(runs(i).time, ' ');
    runs(i).time = runs(i).time((temp.spaces(2)+1):end);
    runs(i).time = datenum( runs(i).time );
    %did this run do beam interpolation?:
    runs(i).bi   = xlsread( [folderstring fileStruct(i).name], 'CTRL', 'B46' );
    runs(i).standard = fileStruct(i).name(1:6); % NBS981 or NBS982?


    %% pick out the '20X' bit of the strings, match to ratios or intensities
    parsei.isosHere = zeros(temp.numfun,2); % mass Table, with intensities in col 1, ratios as (col 1)/(col 2)
    for j = 1:temp.numfun
        parsei.namej = parsei.nameCell{j};
        parsei.start20 = strfind(parsei.nameCell{j}, '20'); %start locations of isotope information
        for k = 1:size(parsei.start20,2)
            parsei.rangek = parsei.start20(k):(parsei.start20(k)+2);
            parsei.isox = str2double(parsei.namej(parsei.rangek)); %pull out each isotope 
            parsei.incell = find(~(parsei.isosHere(j,:) ~= parsei.isox)); %position of isotope in array, if here
            if isempty(parsei.incell) %if the isotope does not already appear
                %place isotope in new column:
                parsei.newcolnum =  nnz(parsei.isosHere(j,:)) + 1; %number of non-zero elements plus one
                parsei.isosHere(j,parsei.newcolnum) = parsei.isox;
            end %-                     if isotope does not appear
        end %for k = instances of an isotope in the name of an output
    end 
    
    %% find intensity measurements, alert if more than one potential
    temp.foundv = zeros(1, temp.numfun); % e.g. [i204 i205 i206 i207 i208]
    for j = 1:temp.numfun
        temp.rowj = parsei.isosHere(j,:);
        if temp.rowj(1)==204 && temp.rowj(2)==0
            if temp.foundv(1) == 1
                disp('already found an i204')
            end
            index.i204 = j;
            temp.foundv(1) = 1;
        elseif temp.rowj(1)==205 && temp.rowj(2)==0
            if temp.foundv(2)==1
                disp('already found an i205')
            end
            index.i205 = j;
            temp.foundv(2) = 1;
        elseif temp.rowj(1)==206 && temp.rowj(2)==0
            if temp.foundv(3)==1
                disp('already found an i206')
            end
            index.i206 = j;
            temp.foundv(3) = 1;
        elseif temp.rowj(1)==207 && temp.rowj(2)==0
            if temp.foundv(4)==1
                disp('already found an i207')
            end
            index.i207 = j;
            temp.foundv(4) = 1;
        elseif temp.rowj(1)==208 && temp.rowj(2)==0
            if temp.foundv(5)==1
                disp('already found an i208')
            end
            index.i208 = j;
            temp.foundv(5) = 1;
        end %if rows
    end %for all rows of isosHere
     
    
    %% now go get the data
    
    % dimensions of cycle data
    temp.ratiosSoFar = xlsread( [folderstring fileStruct(i).name], 'CTRL', 'D17' ); %total # of cycles
    temp.startRowData = sum(isnan(temp.numdata(:,1))) + 1; %there is an extra row when a variable name is numeric.
    temp.cycleRows = temp.startRowData:(temp.startRowData + temp.ratiosSoFar - 1);
    
    % raw data for     [i204   i205   i206   i207   i208]
    runs(i).dataRaw = [temp.numdata(temp.cycleRows, index.i204+2) ...
                       temp.numdata(temp.cycleRows, index.i205+2) ...
                       temp.numdata(temp.cycleRows, index.i206+2) ...
                       temp.numdata(temp.cycleRows, index.i207+2) ...
                       temp.numdata(temp.cycleRows, index.i208+2)]; %plus 2 to make room for cycle# and time
     
	%undo dead time correction, so that dt = 0 for dataDT0 
    runs(i).dataDT0 = runs(i).dataRaw ./ (1 + runs(i).dt*10^-9 * runs(i).dataRaw);
    
    
    %% finally, record cycle times and (linearly interpolated) temperatures
    % where 'secs' are the seconds elapsed since the analysis began and
    % temp is the linearly interpolated temperature from MONITORING worksheet
    
    runs(i).secs = temp.numdata(temp.cycleRows, 2); %time stamp for each cycle
    temp.cyclenums = temp.numdata(temp.cycleRows, 1); 
    temp.monitoring = xlsread( [folderstring fileStruct(i).name], 'MONITORING');
    % temperature data for [cycle num, pyrometer temp]
    temp.tempdata = temp.monitoring(:,[1 7]);
    if temp.ratiosSoFar > temp.tempdata(end,1)
        temp.tempdata(end+1,:) = [temp.ratiosSoFar temp.tempdata(end, 2)];
    end
    runs(i).temp = interp1(temp.tempdata(:,1), temp.tempdata(:,2), 1:temp.ratiosSoFar)';
    
end

runs = runs';  %make each run a row