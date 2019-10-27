function [xref, sref] = cSMCmarg(y,N,xref, prior, sxref, syref,as)
% Generate new refernce trajectory using marginalized cSMC (keep conditioned path on last place) 
% Input: model parameters (prior), observations (y), reference trajectory (xref), nr particles (N)
%, statistics for reference trajectory (sxref,syref), ancestor sampling (as, do if as=1)
% output: new reference trajectory (xref) and its statistics (sref)

% PARAMETERS
%%
[~,T] = size(y);
threshold = 1.0; % For ESS

% Prior hyperparameters
a0 = prior(1);
b0 = prior(2)*ones(N,1);
c0 = prior(3);
d0 = prior(4)*ones(N,1);


% DEFINE VARIABLES
%%
x = zeros(N,T);
w = zeros(N,T);
wun = zeros(N,T);
anc = zeros(N,T);
sx = zeros(N,1);
sy = zeros(N,1);


% INITIALIZE
%%

x0 = zeros(N,1); % Initial states
x0(N) = xref(1); % cSMC
w0 = log(ones(N,1)/N); % initial weights (log)

% Hyperarameters 
alphax = a0;
betax = b0;
alphay = c0;
betay = d0;

% Hyperparameters reference trajectory (for AS)
alphaxT = a0 + T/2;
betaxT = zeros(N,1);
alphayT = c0 + T/2;
betayT = zeros(N,1);


% t=1
%%

% RESAMPLE
Ess = round(ESS(w0));   
if Ess <= threshold*N
    if as==1
        % Update statistics reference trajectory
        sxref = sxref - s(xref(2),f(xref(1),1));
        syref = syref - s(y(1),o(xref(2)));
        
        anc(:,1) = resample(w0);
        
        % Update hyperparameters reference trajectory
        betaxT = betax -s(xref(2), f(x0,1)) - sxref;
        betayT = betay -s(y(1),o(xref(2))) - syref;
        
        wASun = w0 + g(alphax,betax) + g(alphay,betay) - g(alphaxT,betaxT) - g(alphayT,betayT);
        
        % Fix (if NAN weights)
        wASun(betay == Inf) = - Inf;

        wAS = normalize(wASun);
        anc(N,1) = sampleAS(wAS);
        wun(:,1) = log(ones(N,1)/N);
        
    else %if as==0
        anc(:,1) = resample(w0);
        anc(N,1) = N; 
        wun(:,1) = log(ones(N,1)/N);
    end
else % no resample
    anc(:,1) = 1:N;
    wun(:,1) = w0;
end

% PROPAGATE
x(:,1) = propagateMarg(1,alphax,betax(anc(:,1)),x0(anc(:,1)));
x(N,1) = xref(2);


% WEIGHT+NORMALIZE (log)
wun(:,1) = wun(:,1) + weightMarg(alphay,betay(anc(:,1)),y(1),x(:,1));  
w(:,1) = normalize(wun(:,1));


% UPDATE STATISTICS + HYPERPARAMETERS 
sx = sx(anc(:,1)) + s(x(:,1),f(x0(anc(:,1)),1));
sy = sy(anc(:,1)) + s(y(1),o(x(:,1)));

alphax = alphax + 1/2;
betax = b0 - sx; 
alphay = alphay + 1/2;
betay = d0 - sy;



% t=2:T
%%
for t=2:T
    
    % RESAMPLE
    Ess = round(ESS(w(:,t-1)));
    if Ess <= threshold*N
        
        if as == 1  % ancestor sampling
            
            % Update statistics reference trajectory
            sxref = sxref - s(xref(t+1),f(xref(t),t));
            syref = syref - s(y(t),o(xref(t+1))); 
            
            anc(:,t) = resample(w(:,t-1));
            
            % Update hyperparameters reference trajectory
            betaxT = betax -s(xref(t+1), f(x(:,t-1),t)) - sxref;
            betayT = betay -s(y(t),o(xref(t+1))) - syref;
            
            wASun = w(:,t-1) + g(alphax,betax) + g(alphay,betay) - g(alphaxT,betaxT) - g(alphayT,betayT);
            
            % Fix (if NAN weights)
            wASun(betay == Inf) = - Inf;
            
            wAS = normalize(wASun);
            anc(N,t) = sampleAS(wAS);
            wun(:,t) = log(ones(N,1)/N);
            
        else %if as ==0 % No ancestor sampling
            anc(:,t) = resample(w(:,t-1));
            anc(N,t) = N; 
            wun(:,t) = log(ones(N,1)/N);
        end
    else % No resample
        anc(:,t) = 1:N;
        wun(:,t) = w(:,t-1);
    end
    
    % PROPAGATE
    x(:,t) = propagateMarg(t,alphax,betax(anc(:,t)),x(anc(:,t),t-1));
    x(N,t) = xref(t+1); 
    
    % WEIGHT+NORMALIZE+LIKELIHOOD (log)
    wun(:,t) = wun(:,t) + weightMarg(alphay,betay(anc(:,t)),y(t),x(:,t));
    w(:,t) = normalize(wun(:,t));
    
    % UPDATE SUFFICIENT STATISTICS AND HYPERPARAMETERS PARTICLES
    sx = sx(anc(:,t)) + s( x(:,t),f(x(anc(:,t),t-1),t) );
    sy = sy(anc(:,t)) + s( y(t),o(x(:,t)) );
    alphax = alphax + 1/2;
    betax = b0 - sx; 
    alphay = alphay + 1/2;
    betay = d0 - sy;  
    
end


% Draw b according to w_T and extract trajectory for sample b
[xref,b] = sampleTrajectory(x,w,anc,x0);
sref = [sx(b) sy(b)]; 


end


