function w = weightMarg(alpha,beta,y,x)
% Computes marginalized log weights
% Input: states (x), observations (y), parameter (R)
% Output: log-weights (w)
% Model: y = o(x) + e, e~N(0,R)

mu = o(x); % mean value
nu = 2*alpha; % degrees of freedom
sig2 = beta./alpha; % scale parameter

w = stuTpdflog(y,mu,sig2,nu); % log student T pdf



end