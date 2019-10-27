function xt = propagate(t,Q,xt1)
% Propagates states to the next time step
% Input: next time (t), parameters (Q), current states x_t-1 (xt1)
% Output: states x_t at net time step (xt)

N = length(xt1);
xt = f(xt1,t) + sqrt(Q)*randn(N,1);

end