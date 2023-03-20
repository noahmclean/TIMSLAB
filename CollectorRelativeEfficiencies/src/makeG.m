function G = makeG(d, mhat, s2, m0, tails, setup, B, method, dhat)
%MAKEG Make the linear model matrix
%   derivative of dhat with respect to mhat

nData = length(dhat);
nVars = length(mhat);

G = zeros(nData, nVars);

for iVar = 1:nVars

    mPrime = mhat;
    deltaim = setup.relativeStepSize*mhat(iVar);
    imPrime = mhat(iVar) + deltaim;
    mPrime(iVar) = imPrime;

    dPrime = evaluateModel(d, mPrime, s2, m0, tails, setup, B, method);

    G(:,iVar) = (dPrime - dhat)/deltaim;

end % for iVar

end

