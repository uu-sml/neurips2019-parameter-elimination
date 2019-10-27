function [ess] = ESS(w)
% Computes the ESS
% Input: log weights w
% Output: Effective sample size (ESS)

ess = 1/sum(exp(w).^2);

end

