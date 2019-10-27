function [s] = s(x,mu)
% Computes the statistics for a normal l.hood with unknown variance
s = -0.5*(x-mu).^2;
end

