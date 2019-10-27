function [xtraj,b] = sampleTrajectory(x,w, anc,x0)
% Sample new reference trajectory
% Input: States all particles (x), normalized weights (w), ancestors all
% particles (anc) and initial state (x0)
% Output: New refernce trajectory (xtraj) and it's index (b)
[~,T] = size(x);
u = rand(1);
b = find(u<=cumsum(exp(w(:,T))),1);

xtraj(T+1,1) = x(b,T);
a = anc(b,T);
for t=T:-1:2
    xtraj(t,1) = x(a,t-1);
    a = anc(a,t-1);
end
xtraj(1,1) = x0(a);
end