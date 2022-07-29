function G = makeG(dhat, mhat)
%MAKEG Make the linear model matrix
%   derivative of dhat with respect to mhat

nData = length(dhat);
nVars = length(mhat);

G = zeros(nData, nVars);

end

