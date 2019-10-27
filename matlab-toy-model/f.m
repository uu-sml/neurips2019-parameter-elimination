function xt = f(xt1,t)
% Transition model (no noise)
% Input: state at time t-1 (xt1) and time (t)
% Output: state at new time t (xt)
 xt = 0.5*xt1 + 25*xt1./(1+xt1.^2) + 8*cos(1.2*t);