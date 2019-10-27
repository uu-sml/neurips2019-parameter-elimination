function [g] = g(a,b)
% Computes g, the log normalizing factor och the conjugate prior
% Input: Parameters a,b (for inv gamma)
% Output: log g

g = a*log(b) - log(gamma(a));
end

