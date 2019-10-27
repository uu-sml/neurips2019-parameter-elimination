function [y,x] = GenDataNL(Q,R,T,x0)
% Generates states and observations according to the non-linear gaussian model:
% x(t) = f(x(t-1)) + v(t), v(t)~N(0,Q)
% y(t) = o(x(t)) + w(t), w(t)~N(0,R)
% Input: True parameters (Q,R), number of timesteps (T), initial value (x0)
% Output: Observations (y) and true states (x)


x = zeros(1,T);
y = zeros(1,T);

for t=1:T
    if t==1
        x(:,t) = f(x0,t) + sqrt(Q)*randn(1);      
    else
        x(:,t) = f(x(:,t-1),t) + sqrt(Q)*randn(1);  
    end
    
    y(:,t) = o(x(:,t)) + sqrt(R)*randn(1);   
end


end




