clearvars;
close all;

rng(2);

% GENERATE DATA
%%

% Model
% x(t) = f(x(t-1,t) + v(t)   ,    v(t)~N(0,Q)
% y(t) = o(x(t)) + w(t)      ,    w(t)~N(0,R)


% Variance parameters (true)
Q0 = 10; % variance processnoise
R0 = 1; % variance measurement noise

% Initial condition true model
x0 = 0; 

% Generate data (observations y_1:T)
T = 150;
[y,x] = GenDataNL(Q0,R0,T,x0);         

pk = 500; % Print iteration nr every pk iteration


% SET PARAMETERS MCMC
%%

% Set parameters MCMC (N proportional to T)
Ninit = 2000; % Number of particles initital run (to generate xref)
Npgas = 50; % Number of particles in PGAS (unmarginalized)
Npgasm = 50; % Number of particles in PGAS (marginalized)

M = 10000; % Number of iterations
bi = 150; % Burn in

% Hyperparameters invGamma prior
a = .001;
b = .001;
c = 1;
d = 1;


% INITIALIZE MCMC CHAIN
%%

% Initialize parameters
Qinit = 10^2;
Rinit = 10^2;

% Initial ref trajectory (same all methods)
xrefinit = pf(Qinit,Rinit,y,Ninit);

% Compute statistic initial xref
sxref = 0;
syref = 0;
for k = 1:T
   sxref = sxref + s(xrefinit(k+1),f(xrefinit(k),k)); 
   syref = syref + s(y(k),o(xrefinit(k+1))); 
end
srefinit = [sxref syref];



% SETTINGS BLOCKING
%%
B = 5;
L = 20;


% INFER STATES AND PARAMETERS (NOT MARGINALIZED, PGAS)
%%

as = 1; % Ancestor sampling

tic; % Timer

Qpgas = zeros(1,M);
Rpgas = zeros(1,M);
xrefPGAS = zeros(T+1,M);
update_pgas = zeros(T+1,1);

% Initialize parameters
Qpgas(1) = Qinit;
Rpgas(1) = Rinit;

xrefPGAS(:,1) = xrefinit;

for m=2:M    
    
    % Sample x
    xrefPGAS(:,m) = cSMC(Qpgas(m-1),Rpgas(m-1),y,Npgas,xrefPGAS(:,m-1),as); 
     % Sample theta
    [Qpgas(m), Rpgas(m)] = sampleParam(xrefPGAS(:,m),y,a,b,c,d);
    
    update_pgas(:,1) = update_pgas(:,1) + (xrefPGAS(:,m)~=xrefPGAS(:,m-1)); % Update frequency
    
    if mod(m,pk)==0 % Print iteration nr
        m
    end
end

% Update frequency
update_pgas = update_pgas/M;
 
tpgas = toc; 





% INFER STATES AND PARAMETERS (MARGINALIZED, PGAS)
%%

as=1; % ancestor sampling

tic; % Timer

Qpgasm = zeros(1,M);
Rpgasm = zeros(1,M);
xrefPGASm = zeros(T+1,M);
update_pgasm = zeros(T+1,1);

% Initialize parameters
Qpgasm(1) = Qinit;
Rpgasm(1) = Rinit;

xrefPGASm(:,1) = xrefinit;
srefPGASm = srefinit;



for m=2:M    
    
    % Sample x
    [xrefPGASm(:,m), srefPGASm] = cSMCmarg(y,Npgasm,xrefPGASm(:,m-1),[a;b;c;d], srefPGASm(1),srefPGASm(2),as);
    % Sample theta
    Qpgasm(m) = invGamma(a+T/2, b - srefPGASm(1));
    Rpgasm(m) = invGamma(c+T/2, d - srefPGASm(2));

    
    update_pgasm(:,1) = update_pgasm(:,1) + (xrefPGASm(:,m)~=xrefPGASm(:,m-1));
    
    if mod(m,pk)==0  % Print iteration nr
        m       
    end
end

% Update frequency reference trajectory
update_pgasm = update_pgasm/M;

tpgasm = toc;




% INFER STATES AND PARAMETERS (MARGINALIZED, BLOCKED, PGAS)
%%

as=1; % ancestor sampling
tic;  % Timer

Qpgasmb = zeros(1,M);
Rpgasmb = zeros(1,M);
xrefPGASmb = zeros(T+1,M);
update_pgasmb = zeros(T+1,1);

% Initialize parameters
Qpgasmb(1) = Qinit;
Rpgasmb(1) = Rinit;

xrefPGASmb(:,1) = xrefinit;
srefPGASm = srefinit;



for m=2:M    
    
    xrefPGASmb(:,m) = xrefPGASmb(:,m-1); % copy xref previous step to new step
    
    % Block 1, update ref traj [0:B+L] cond on y[1:B+L] and ref traj [B+L+1]
    xrefPGASmb(1:B+L+1,m) = cSMCblock(Qpgasmb(m-1),Rpgasmb(m-1),y(1:B+L),Npgasm, xrefPGASmb(1:B+L+1,m), xrefPGASmb(B+L+2,m), as);
    
    % Block 2
    [xrefPGASmb(:,m), srefPGASm] = cSMCblockmarg(y,Npgasm,xrefPGASmb(:,m),[a;b;c;d],as,B); 
    
    % Sample theta
    Qpgasmb(m) = invGamma(a+T/2, b - srefPGASm(1));
    Rpgasmb(m) = invGamma(c+T/2, d - srefPGASm(2));
    

    update_pgasmb(:,1) = update_pgasmb(:,1) + (xrefPGASmb(:,m)~=xrefPGASmb(:,m-1));
    
    if mod(m,pk)==0 % Print iteration nr
        m       
    end
end

% Update frequency reference trajectory
update_pgasmb = update_pgasmb/M;


tpgasmb = toc;



 
figure(1);
p = plot(1:T,update_pgas(2:end), 'b', 1:T,update_pgasm(2:end), 'b:', 1:T,update_pgasmb(2:end),'b--');
set(p,'LineWidth',1.5);
xlabel('time');
ylabel('update freq x');
legend('PG','PGm','PGmb','PGAS','PGASm', 'PGASmb');
