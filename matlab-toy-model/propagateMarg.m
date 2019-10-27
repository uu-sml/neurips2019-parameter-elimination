function xt = propagateMarg(t,alpha,beta,xt1)
% Propagates states to the next time step, marginalized
% Input: next time (t), parameters (alpha,beta), current states x_t-1 (xt1)
% Output: states x_t at net time step (xt)

N = length(xt1);
mu = f(xt1,t); % mean value
nu = 2*alpha; % degrees of freedom
sig2 = beta./alpha; % scale parameter

xt = mu + sqrt(sig2).*trnd(nu,N,1);

end