function s2 = calculateUnct(data, d, method, setup)
%CALCULATEUNCT a priori measurement uncertainty
%   Calculate using Johnson noise and shot noise estimates
%   Take integration periods from parsed method file
%   Use amplifier data from parssed data file header/COLLECTORS block
%   Noah McLean for faradayRelativeEfficiency, started on 25 July 2022

%% grab integration period from each baseline and on-peak sequence

nBLseq = size(method.baselines,2);
nOPseq = size(method.onpeaks,2);

integPeriodBL = zeros(nBLseq,1);
for iBL = 1:nBLseq

    integPeriodString = string(method.baselines(iBL).Info(5).Value); % e.g. "ms100"
    integPeriodms = extractAfter(integPeriodString, digitBoundary); % extract milliseconds as number
    integPeriodBL(iBL) = double(integPeriodms);

end % for iBL

integPeriodOP = zeros(nOPseq,1);
for iOP = 1:nOPseq

    integPeriodString = string(method.onpeaks(iOP).Info(6).Value); % e.g. "ms100"
    integPeriodms = extractAfter(integPeriodString, digitBoundary); % extract milliseconds as number
    integPeriodOP(iOP) = double(integPeriodms);

end % for iBL


%% make a vector of length d that contains integration periods
% integPeriod is integration period of each element in d, in seconds

isBL = ~d.isOP;
integPeriod = zeros(size(d.int));
BLseqs = d.seq .* isBL;
OPseqs = d.seq .* d.isOP;
integPeriod(isBL)   = integPeriodBL(BLseqs(BLseqs>0)) * 0.001; % in seconds
integPeriod(d.isOP) = integPeriodOP(OPseqs(OPseqs>0)) * 0.001; 


%% calculate Johnson noise

% Johnson noise variance = 4 * kB * T * R * frequency
dResist = data.FaradayResist(d.det); % resistance, same size as d
s2 = 4 * setup.kB * setup.tempInK * dResist ./ integPeriod;


%% calculate shot noise to On Peak 

shotNoise = d.int .* (dResist/setup.coulomb) ./ integPeriod;
s2 = s2 + shotNoise .* d.isOP;

end % function caculateUnct
