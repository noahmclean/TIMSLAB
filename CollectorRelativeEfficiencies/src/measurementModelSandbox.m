n = 1e6;
a = random('normal', 1, 0.01, [n 1]);
b = random('normal', 2, 0.014, [n 1]);
s2 = [0.01^2*ones(n,1); 0.014^2*ones(n, 1)];

alrab = log(a./b);
logMean = mean(alrab);
logStde = std(alrab)/sqrt(length(alrab));

d = [a; b];

m0 = [log(1/2); log(2)];

% check that model parameterization works
dhat = evalg(a, b, m0);

% best fit to data using fminunc
opts = optimoptions("fminunc");
opts.StepTolerance = 1e-10;
[x, fval] = fminunc(@(m) objfunc(d, a, b, s2, m), m0, opts);

%calculate uncertainty
G = makeG(a, b, x);
covx = inv( (G.*(1./s2))' * G ); % inv(G'*W*G), where W = diag(1./s2);
unctx = sqrt(diag(covx));

% x(1) = alrab -- the best fit is the alr mean!
% unctx(1) == logStde when all(a > 0 & b > 0).  Awesome!

% Note after coding: this all works out!  No logs of negative ratios
% due to negative measured intensities (a and b)!! Best fit value 
% of model parameters converges to geometric mean/alr results
% and uncertainty converges to standard error of logratios
% when all(a > 0 & b > 0)



%% functions

function G = makeG(a, b, m)

    arange = 1:length(a);
    brange = max(arange)+1 : max(arange) + length(b);
    
    lograbtrue = m(1);
    logbtrue = m(2);

    G(arange,1) = exp(lograbtrue + logbtrue);
    G(arange,2) = exp(lograbtrue + logbtrue);
    G(brange,1) = zeros(size(b));
    G(brange,2) = exp(logbtrue);

end % fuction makeG


% evaluate model function g(m), 
% return predicted data vector dhat
function dhat = evalg(a, b, m)

    arange = 1:length(a);
    brange = max(arange)+1 : max(arange) + length(b);
    
    dhat = zeros(size([a; b]));
    dhat(arange) = exp(m(1) + m(2));
    dhat(brange) = exp(m(2));

end % function evalg


% objective function - weighted sum of squares for d-dhat
function fval = objfunc(d, a, b, s2, m)

    dhat = evalg(a, b, m);
    r = d - dhat;
    fval = sum( r.^2 ./ s2);

end %function objfunc
