function chi2 = objfunc(d, m, s2, m0, tails, setup, B, method)
%OBJFUNC Objective function for faraday efficiency fit
%   return chi square statistic (not reduced)


dhat = evaluateModel(d, m, s2, m0, tails, setup, B, method);

r = d.int - dhat;
rejects = abs(r) > 0.5;

% remove major outliers (needs investigation)
rClean = r(~rejects);
s2Clean = s2(~rejects);

chi2 = sum( rClean.^2 ./ s2Clean );


end % function objfunc

