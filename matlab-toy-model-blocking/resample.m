function anc = resample(w)
% Multiomial resampling
% Input: log-weights (w)
% Output: N ancestor indicies (anc) for surviving particles

N = length(w);
c = cumsum(exp(w));
anc = zeros(N,1);
for i=1:N
    u = rand(1);
    anc(i) = find(u<=c,1);
end

end