function [xref] = cSMCblock(Q,R,y,N,xref, xfinal, as)
% Generate new refernce trajectory for block 1 using cSMC (keep conditioned path on last place)
% Input: model parameters (Q,R), observations (y), reference
% trajectory (xref), nr particles (N), last state in block 1 (xfinal), ancestor sampling (as, do if as=1)
% output: new reference trajectory (xref)

% PARAMETERS
%%
[~,T] = size(y);
threshold = 1;


% INITIALIZE  
%%
x0 = zeros(N,1); % Initial states
x0(N) = xref(1); % cSMC
w0 = log(ones(N,1)/N); % initial weights (log)

x = zeros(N,T);
w = zeros(N,T);
wun = zeros(N,T);
anc = zeros(N,T);


% FOR t=1
%%
% RESAMPLE
Ess = ESS(w0);
if Ess <= threshold*N
    anc(:,1) = resample(w0);
    if as==1
        wun_as = w0 - 0.5*log(2*pi*Q) - 0.5*(xref(2)-f(x0,1)).^2./Q;
        wAS = normalize(wun_as);
        anc(N,1) = sampleAS(wAS); 
    else
        anc(N,1) = N; 
    end
    wun(:,1) = ones(N,1)/N;
else
    anc(:,1) = 1:N;
    wun(:,1) = w0;
end

% PROPAGATE
x(:,1) = propagate(1,Q,x0); 
x(N,1) = xref(2); 


% WEIGHT+NORMALIZE (log)
wun(:,1) = weight(R,y(1),x(:,1));        
w(:,1) = normalize(wun(:,1));


% FOR t=2:T
%%
for t=2:T
    
    % RESAMPLE
    Ess = ESS(w(:,t-1));
    if Ess <= threshold*N
        anc(:,t) = resample(w(:,t-1));
        if as==1
            wun_as = w(:,t-1) - 0.5*log(2*pi*Q) - 0.5*(xref(t+1)-f(x(:,t-1),t)).^2./Q; 
            wAS = normalize(wun_as);
            anc(N,t) = sampleAS(wAS); 
        else
            anc(N,t) = N; 
        end
        wun(:,t) = log(ones(N,1)/N);
    else
        anc(:,t) = 1:N;
        wun(:,t) = w(:,t-1);
    end
    
    % PROPAGATE
    x(:,t) = propagate(t,Q,x(anc(:,t),t-1)); 
    x(N,t) = xref(t+1); 
    
    % WEIGHT+NORMALIZE+LIKELIHOOD (log)
    wun(:,t) = wun(:,t) + weight(R,y(t),x(:,t)); 
    if t==T
        wun(:,t) = wun(:,t) + weightFinal(Q,xfinal,x(:,t),t+1);  
    end
    w(:,t) = normalize(wun(:,t));
       
end


% FINISH
% Draw b according to w_T and extract trajectory for sample b
[xref,b] = sampleTrajectory(x,w,anc,x0);



end


