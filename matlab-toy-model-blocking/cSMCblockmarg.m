function [xref, sref] = cSMCblockmarg(y,N,xref, prior, as,B)
% Generate new refernce trajectory for block 2 using marginalized cSMC (keep conditioned path on last place) 
% Input: model parameters (prior), observations (y), reference trajectory
% (xref), nr particles (N), ancestor sampling (as, do if as=1)
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

x0 = xref(1)*ones(N,1); % Initial states (from ref traj)

x = ones(N,1)*xref(2:end)'; % Keep reference trajectory for first B timesteps


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

% Statistic complete reference trajectory
sxref = 0;
syref = 0;
for k = 1:T
    sxref = sxref + s(xref(k+1),f(xref(k),k));
    syref = syref + s(y(k),o(xref(k+1)));
end


% for t=1:T
%%

for t=1:T
    
    if t<=B % Only compute statistics up to time B
        
        if as == 1  % ancestor sampling   
            % Update statistics reference trajectory
            sxref = sxref - s(xref(t+1),f(xref(t),t)); 
            syref = syref - s(y(t),o(xref(t+1))); 
            
            % Update hyperparameters reference trajectory
            if t==1
                betaxT = betax -s(xref(t+1), f(xref(t),t)) - sxref;
            else
                betaxT = betax -s(xref(t+1), f(x(:,t-1),t)) - sxref;
            end
            betayT = betay -s(y(t),o(xref(t+1))) - syref;
        end

        % No resample
        anc(:,t) = 1:N;
        wun(:,t) = log(ones(N,1)/N);
        w(:,t) = normalize(wun(:,t));
        
        % Update statistics + hyperparameters
        if t==1
            sx = sx(anc(:,t)) + s( x(:,t),f(xref(t),t) );
        else
            sx = sx(anc(:,t)) + s( x(:,t),f(x(anc(:,t),t-1),t) );
        end
        sy = sy(anc(:,t)) + s( y(t),o(x(:,t)) );
        alphax = alphax + 1/2;
        betax = b0 - sx;
        alphay = alphay + 1/2;
        betay = d0 - sy;
        
    else % For t>B, run marginalized cSMC
        % RESAMPLE
        Ess = round(ESS(w(:,t-1)));
        if Ess <= threshold*N
            
            if as == 1  % ancestor sampling
                
                % UUpdate statistics reference trajectory
                sxref = sxref - s(xref(t+1),f(xref(t),t)); 
                syref = syref - s(y(t),o(xref(t+1))); 
                
                anc(:,t) = resample(w(:,t-1));
                
                % Update hyperparameters reference trajectory
                betaxT = betax -s(xref(t+1), f(x(:,t-1),t)) - sxref;
                betayT = betay -s(y(t),o(xref(t+1))) - syref;
                
                wASun = w(:,t-1) + g(alphax,betax) + g(alphay,betay) - g(alphaxT,betaxT) - g(alphayT,betayT);
                
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
        
        % WEIGHT+NORMALIZE (log)
        wun(:,t) = wun(:,t) + weightMarg(alphay,betay(anc(:,t)),y(t),x(:,t));
        w(:,t) = normalize(wun(:,t));
        
        % UPDATE SUFFICIENT STATISTICS AND HYPERPARAMETERS
        sx = sx(anc(:,t)) + s( x(:,t),f(x(anc(:,t),t-1),t) );
        sy = sy(anc(:,t)) + s( y(t),o(x(:,t)) );
        alphax = alphax + 1/2;
        betax = b0 - sx;
        alphay = alphay + 1/2;
        betay = d0 - sy;
     
    end
end
 




% Draw b according to w_T and extract trajectory for sample b
[xref,b] = sampleTrajectory(x,w,anc,x0);
sref = [sx(b) sy(b)]; 


end


