function p = invGamma(a,b)
% Function to generate a sample from the inverse gamma distribution
% Input: Parameters a (shape) and b (scale)
% Output: sample from distribution

p = 1./gamrnd(a,1./b);

end