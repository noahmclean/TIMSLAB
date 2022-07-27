function chi2 = objfunc(d, m, s2, m0, tails, setup)
%OBJFUNC Objective function for faraday efficiency fit
%   return chi square statistic (not reduced)


dhat = evaluateModel(d, m, m0, tails, setup);

chi2 = sum( (d.int - dhat)./s2 );


end % function objfunc
