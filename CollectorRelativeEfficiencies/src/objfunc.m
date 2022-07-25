function chi2 = objfunc(d, s2, m, m0)
%OBJFUNC Objective function for faraday efficiency fit
%   return chi square statistic (not reduced)


dhat = evaluateModel(d, m, m0);

chi2 = sum( (d.int - dhat)./s2 );


end % function objfunc

