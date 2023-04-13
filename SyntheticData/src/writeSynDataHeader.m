function writeSynDataHeader(setup, filename)
%WRITEHEADER Write the header for a synthetic data file
%   Copy header format from Version 2.0.11 of Isolynx, 2022
%
% Written for Tripoli and CollectorRelativeEfficiencies on 5-Apr-2023
% by Noah McLean



% if  writeData is false, return to makeSyntheticData
if ~setup.writeData
    return
end

% create header block
methodName = setup.method.methodName + ".TIMSAM";

header = ...
    ["#HEADER";
     "Analysis";
     "Version,"        + "2.0.11,1.04"; 
     "Filename,"       + filename + ".TXT";
     "MethodName,"     + methodName;
     "MethodPath,"     + "/TIMSLAB/SyntheticData/Methods";
     "IsoWorksMethod," + "This is synthetic data from the TIMSLAB repository";
     "FolderPath,"     + "/TIMSLAB/SyntheticData/syndata";
     "Corrected,"      + "Yes"; % Presumably corrected for collector gain and efficiency
     "BChannels,"      + "No"; % no ATONA BChannels yet
     "TimeZero,"       + setup.timeCreated;
     "";
     "#COLLECTORS";
     "Name,Type,Resistor,Gain,Efficiency,DT"];

writematrix(header, "../syndata/"+filename, "QuoteStrings", "none")


end % function

