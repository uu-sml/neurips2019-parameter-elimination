function [Q,R] = sampleParam(xref,y,a,b,c,d)
% Samples new parameters from posterior p(theta|x,y)

T = length(xref)-1;
diff_Q = zeros(1,T);
diff_R = zeros(1,T);
for k = 1:T
    diff_Q(k) = xref(k+1) - f(xref(k),k);
    diff_R(k) = y(k)-o(xref(k+1));
end
Q = invGamma(a+T/2, b + 0.5*sum((diff_Q).^2));
R = invGamma(c+T/2, d + 0.5*sum((diff_R).^2));

end

