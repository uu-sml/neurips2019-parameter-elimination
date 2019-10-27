function w = weight(R,y,x)
% Computes log weights
% Input: states (x), observations (y), parameter (R)
% Output: log-weights (w)
% Model: y = o(x) + e, e~N(0,R)

w = -0.5*log(2*pi*R) - 0.5*(y-o(x)).^2./R; % log normal pdf

end