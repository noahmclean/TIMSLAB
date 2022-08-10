%% Reduce data from multidynamic Pb method for measuring Pb peak tails
% PbPeakTails.TIMSAM on Phoenix
%  Needs a reference voltage measurement (LOS closed electronic baselines)
%  before and after the tails measurement (on a 208Pb spike, preferably)

refVolts = mean(refVolts1(450:500,:));
refVolts1sAbs = std(refVolts1(450:500,:));

faraCorr = fara - refVolts;
isSeq1 = seq == "OP1";
isSeq2 = seq == "OP2";
isSeq3 = seq == "OP3";
isSeq4 = seq == "OP4";

r204208S4 = faraCorr(isSeq4,1)./faraCorr(isSeq4,5);
r205208S4 = faraCorr(isSeq4,2)./faraCorr(isSeq4,5);
r206208S4 = faraCorr(isSeq4,3)./faraCorr(isSeq4,5);
r207208S4 = faraCorr(isSeq4,4)./faraCorr(isSeq4,5);
r2095208S4 = faraCorr(isSeq4,6)./faraCorr(isSeq4,5);

r2065208S3 = faraCorr(isSeq3,5)./faraCorr(isSeq3,6);
r2075206S2 = faraCorr(isSeq2,6)./faraCorr(isSeq2,5);
r2075208 = mean(r2075206S2(600:6000))*mean(r206208S4(300:3000));