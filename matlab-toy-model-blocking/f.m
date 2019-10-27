function xt1 = f(xt,t)
% Transition model (no noise)
% Input: state at time t-1 (xt1) and time (t)
% Output: state at new time t (xt)
 xt1 = 0.5*xt + 25*xt./(1+xt.^2) + 8*cos(1.2*t);