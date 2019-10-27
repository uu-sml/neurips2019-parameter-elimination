function w = weightFinal(Q,xfinal,x,t)
% computes weights final time step
% Input: states x, final state xfinal, parameter Q, time t
% Output: log-weights w
% Model: y = o(x) + e, e~N(0,R)

w = -0.5*log(2*pi*Q) - 0.5*(xfinal-f(x,t)).^2./Q; % log normal pdf

end