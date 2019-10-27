function w = normalize(wun) 
% function to normalize log-weights
% Input: unnormalized log-weights wun
% Output: normalized log-weights w

maxw = max(wun);
w_shift = wun - maxw;
w = w_shift - log(sum(exp(w_shift)));

end