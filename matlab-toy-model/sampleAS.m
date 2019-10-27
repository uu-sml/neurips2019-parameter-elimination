function [anc] = sampleAS(w)
% Sample one ancestor index
u = rand(1);
anc = find(u<=cumsum(exp(w)),1);
end

