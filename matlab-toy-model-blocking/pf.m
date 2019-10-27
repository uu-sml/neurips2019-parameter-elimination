function [xref] = pf(Q,R,y,N)
% Standard (bootstrap) particle filter 
% Input: model parameters (Q,R), observations (y), nr particles (N)
% output: refernce trajectory 

% PARAMETERS
%%
[~,T] = size(y);
threshold = 1;


% INITIALIZE  
%%
x0 = zeros(N,1);
w0 = ones(N,1)/N; 

x = zeros(N,T);
w = zeros(N,T);
wun = zeros(N,T);
anc = zeros(N,T);


% FOR t=1
%%
% RESAMPLE 
anc(:,1)=1:N;
w(:,1)=w0;

% PROPAGATE
x(:,1) = propagate(1,Q,x0); 


% WEIGHT+NORMALIZE+LIKELIHOOD (log)
wun(:,1) = weight(R,y(1),x(:,1));         
w(:,1) = normalize(wun(:,1));


% FOR t=2:T
%%
for t=2:T
    
    % RESAMPLE
    Ess = ESS(w(:,t-1));
    if Ess <= threshold*N
        anc(:,t) = resample(w(:,t-1));
        wun(:,t) = log(ones(N,1)/N);
    else
        anc(:,t)=1:N;
        wun(:,t)=w(:,t-1);
    end
    
    % PROPAGATE
    x(:,t) = propagate(t,Q,x(anc(:,t),t-1)); 
    
    % WEIGHT+NORMALIZE+LIKELIHOOD (log )
    wun(:,t) = wun(:,t) + weight(R,y(t),x(:,t));        
    w(:,t) = normalize(wun(:,t));
        
end


% FINISH
% Draw b according to w_T and extract trajectory for sample b
xref = sampleTrajectory(x,w,anc,x0);


end


