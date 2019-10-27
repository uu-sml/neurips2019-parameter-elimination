clearvars;
close all;

rng(2);

% GENERATE DATA
%%

% Model
% x(t) = f(x(t-1),t) + v(t)   ,    v(t)~N(0,Q)
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
Npg = 500;  % Number of particles in PG (unmarginalized)
Npgm = 500; % Number of particles in PG (marginalized)
Npgas = 500; % Number of particles in PGAS (unmarginalized)
Npgasm = 500; % Number of particles in PGAS (marginalized)

M = 10000; % Number of iterations
bi = 1500; % Burn in
lag = 40; % Lag for autocorrelation plots

% Hyperparameters invGamma prior
a = 1; % shape parameter transition
b = 1; % scale paraeter transition
c = 1; % shape parameter observation
d = 1; % scale parameter observation



% INITIALIZE MCMC CHAIN
%%

% Initialize parameters
Qinit = 10^2;
Rinit = 10^2;

% Initial ref trajectory (same all methods)
xrefinit = pf(Qinit,Rinit,y,Ninit);


% INFER STATES AND PARAMETERS (PG, UNMARGINALIZED)
%%

as = 0; % No ancestor sampling

tic; % Timer

Qpg = zeros(1,M);
Rpg = zeros(1,M);
xrefPG = zeros(T+1,M);
update_pg = zeros(T+1,1);

% Initialize
Qpg(1) = Qinit;
Rpg(1) = Rinit;
xrefPG(:,1) = xrefinit;


for m=2:M    
    
    % Sample x
    xrefPG(:,m) = cSMC(Qpg(m-1),Rpg(m-1),y,Npg,xrefPG(:,m-1),as); 
    % SAmple theta
    [Qpg(m), Rpg(m)] = sampleParam(xrefPG(:,m),y,a,b,c,d);  
    
    
    update_pg(:,1) = update_pg(:,1) + (xrefPG(:,m)~=xrefPG(:,m-1)); % Update frequency
    
    if mod(m,pk)==0 % Print iteration nr
        m 
    end
end

% Update frequency 
update_pg = update_pg/M;

% Autocorrelation
acf_Rpg = acf(Rpg(bi:end)',lag);
acf_Qpg = acf(Qpg(bi:end)',lag);

tpg = toc; % Timer




% INFER STATES AND PARAMETERS (PG, MARGINALIZED)
%%

as=0; % No ancestor sampling

tic; % Timer

Qpgm = zeros(1,M);
Rpgm = zeros(1,M);
xrefPGm = zeros(T+1,M);
update_pgm = zeros(T+1,1);

% Initialize 
Qpgm(1) = Qinit;
Rpgm(1) = Rinit;
xrefPGm(:,1) = xrefinit;

% Compute statistic initial refernce trajectory
sxref = 0;
syref = 0;
for k = 1:T
   sxref = sxref + s(xrefPGm(k+1,1),f(xrefPGm(k,1),k)); 
   syref = syref + s(y(k),o(xrefPGm(k+1,1))); 
end
sref = [sxref syref];

for m=2:M   
    
    % Sample x
    [xrefPGm(:,m), sref] = cSMCmarg(y,Npgm,xrefPGm(:,m-1),[a;b;c;d],sref(1), sref(2),as); 
    % Sample theta
    Qpgm(m) = invGamma(a+T/2, b - sref(1));
    Rpgm(m) = invGamma(c+T/2, d - sref(2));
    
    update_pgm(:,1) = update_pgm(:,1) + (xrefPGm(:,m)~=xrefPGm(:,m-1)); % Update frequency
    
    if mod(m,pk)==0 % Print iteration nr
        m
    end
end

% Update frequency 
update_pgm = update_pgm/M;
         
% Autocorrelation 
acf_Rpgm = acf(Rpgm(bi:end)',lag);
acf_Qpgm = acf(Qpgm(bi:end)',lag);

tpgm = toc;



% INFER STATES AND PARAMETERS (PGAS, NOT MARGINALIZED)
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
    
    if mod(m,pk)==0  % Print iteration nr
        m
    end
end

% Update frequency
update_pgas = update_pgas/M;

% Autocorrelation
acf_Rpgas = acf(Rpgas(bi:end)',lag);
acf_Qpgas = acf(Qpgas(bi:end)',lag);
 
tpgas = toc; 


% INFER STATES AND PARAMETERS (PGAS, MARGINALIZED)
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

% Compute statistic initial xref
sxref = 0;
syref = 0;
for k = 1:T
   sxref = sxref + s(xrefPGASm(k+1,1),f(xrefPGASm(k,1),k)); 
   syref = syref + s(y(k),o(xrefPGASm(k+1,1))); 
end
sref = [sxref syref];



for m=2:M    
 
    % Sample x
    [xrefPGASm(:,m), sref] = cSMCmarg(y,Npgasm,xrefPGASm(:,m-1),[a;b;c;d], sref(1),sref(2),as); 
    % Sample theta
    Qpgasm(m) = invGamma(a+T/2, b - sref(1));
    Rpgasm(m) = invGamma(c+T/2, d - sref(2));

    update_pgasm(:,1) = update_pgasm(:,1) + (xrefPGASm(:,m)~=xrefPGASm(:,m-1)); % Update frequency
    
    if mod(m,pk)==0   % Print iteration nr
        m       
    end
end

% Update frequency 
update_pgasm = update_pgasm/M;

% Autocorrelation
acf_Rpgasm = acf(Rpgasm(bi:end)',lag);
acf_Qpgasm = acf(Qpgasm(bi:end)',lag);

tpgasm = toc;





% PLOT RESULTS
%%

Qlow = 0;
Qhigh = 20;
Rlow = 0;
Rhigh = 3.5;

% Trace plots
figure(1);
title('Trace parameters');
subplot(2,1,1);
plot(1:M, Qpg, 1:M, Qpgas, 1:M, Qpgm, 1:M, Qpgasm);
ylabel('Q');
ylim([Qlow,1.25*Qhigh]);
legend('PG','PGAS','PGm','PGASm');
subplot(2,1,2);
plot(1:M, Rpg, 1:M, Rpgas, 1:M, Rpgm, 1:M, Rpgasm);
ylabel('R');
ylim([Rlow,1.25*Rhigh]);
legend('PG','PGAS','PGm','PGASm');

% Histograms
figure(2);
subplot(2,4,1);
histogram(Qpg(bi:end));
hold on;
yLimit = get(gca,'YLim');
line([Q0, Q0], yLimit, 'LineWidth', .5, 'Color', 'k');
hold off;
xlim([Qlow,Qhigh]);
ylabel('Q');
title('PG');
subplot(2,4,5);
histogram(Rpg(bi:end));
hold on;
yLimit = get(gca,'YLim');
line([R0, R0], yLimit, 'LineWidth', .5, 'Color', 'k');
hold off;
xlim([Rlow,Rhigh]);
ylabel('R');
subplot(2,4,2);
histogram(Qpgas(bi:end));
hold on;
yLimit = get(gca,'YLim');
line([Q0, Q0], yLimit, 'LineWidth', .5, 'Color', 'k');
hold off;
xlim([Qlow,Qhigh]);
ylabel('Q');
title('PGAS');
subplot(2,4,6);
histogram(Rpgas(bi:end));
hold on;
yLimit = get(gca,'YLim');
line([R0, R0], yLimit, 'LineWidth', .5, 'Color', 'k');
hold off;
xlim([Rlow,Rhigh]);
ylabel('R');
subplot(2,4,3);
histogram(Qpgm(bi:end));
hold on;
yLimit = get(gca,'YLim');
line([Q0, Q0], yLimit, 'LineWidth', .5, 'Color', 'k');
hold off;
xlim([Qlow,Qhigh]);
ylabel('Q');
title('PGm');
subplot(2,4,7);
histogram(Rpgm(bi:end));
hold on;
yLimit = get(gca,'YLim');
line([R0, R0], yLimit, 'LineWidth', .5, 'Color', 'k');
hold off;
xlim([Rlow,Rhigh]);
ylabel('R');
subplot(2,4,4);
histogram(Qpgasm(bi:end));
hold on;
yLimit = get(gca,'YLim');
line([Q0, Q0], yLimit, 'LineWidth', .5, 'Color', 'k');
hold off;
xlim([Qlow,Qhigh]);
ylabel('Q');
title('PGASm');
subplot(2,4,8);
histogram(Rpgasm(bi:end));
hold on;
yLimit = get(gca,'YLim');
line([R0, R0], yLimit, 'LineWidth', .5, 'Color', 'k');
hold off;
xlim([Rlow,Rhigh]);
ylabel('R');

% Update frequency for refernce trajectory
figure(3);
plot(0:T,update_pg, 0:T,update_pgas,0:T,update_pgm, 0:T,update_pgasm);
xlabel('time');
ylabel('update freq x');
legend('PG','PGAS','PGm','PGASm');

% Autocorrelation for parameters
figure(4);
title('Autocorrelation parameters');
subplot(2,1,1);
plot(0:lag,acf_Qpg,'r',0:lag,acf_Qpgas, 'b',0:lag,acf_Qpgm, 'r:', 0:lag,acf_Qpgasm,'b:');
xlim([0,lag]);
xlabel('lag');
ylabel('Q');
legend('PG','PGAS','PGm','PGASm');
subplot(2,1,2);
plot(0:lag,acf_Rpg, 0:lag,acf_Rpgas,0:lag,acf_Rpgm, 0:lag,acf_Rpgasm);
xlim([0,lag]);
xlabel('lag');
ylabel('R');
legend('PG','PGAS','PGm','PGASm');

% Running time
tpg
tpgm
tpgas
tpgasm

