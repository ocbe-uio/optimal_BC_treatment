%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to simulate alternating treatment
%
% Author: Tuğba Akman Date: Jan 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data are taken from
% Obesity-Activated Adipose-Derived Stromal Cells Promote
% Breast Cancer Growth and Invasion
% Neoplasia (2018) 20, 1161-1174
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main_AlternatingTreatment()

close all
clear
clc
format long

set(0, 'defaultaxesfontsize',16)

global Eq


%Initial conditions
T0=1;
S0=0.75;
R0=T0-S0;
E0_CD=170;
E0_HFD=1200;

F0_CD= 50;
F0_HFD = 360;

% Parameters
Eq.mu = 5.94;
Eq.m1 = 5e-4;
Eq.m2CD = 1/(F0_HFD);
Eq.m2HFD = 1/(F0_HFD);
Eq.n = 1;
Eq.k1=0.55;
Eq.a1=43;
Eq.r=20;
Eq.alpha=1.7e-6;
Eq.pCD=1;
Eq.pHFD=1;
Eq.K=1/Eq.m1;
Eq.m=1;
Eq.betaCD=1;
Eq.betaHFD=1;
Eq.l=10;
Eq.Sbound=0;
Eq.k2CD = 0.05;
Eq.k2HFD= 0.05;
Eq.k3CD=0.5*Eq.k1;
Eq.k3HFD=0.5*Eq.k1;
Eq.a2=10;
Eq.a3=10;

% Weights in the cost functional
Eq.weight1 = 1;
Eq.weight2 = 1;
Eq.weight3 = 1;
Eq.weight4 = 1;
Eq.weight5 = 1;
Eq.weight6 = 1;

% % Upper bound for control or max treatment
M2=0.99;            

% Time discretization
t0=0; % initial time
tf = 25; % final time
Interval=100; % Number of subintervals
dt = ( tf-t0)/Interval; % step size
Tu=t0:(1*dt):tf;

% tumor_threshold to start treatment
tumor_threshold = ceil(Eq.K*0.25);

initx=[S0 R0 E0_CD F0_CD S0 R0 E0_HFD F0_HFD]';% IC for state

%Allocate for uncontrolled case - MUST BE EQUAL TO ZERO!!!!
uCD=zeros(size(Tu))';
uHFD=zeros(size(Tu))';

%Toilerance for ode solver
options = odeset('RelTol',1e-6,'AbsTol',1e-6);

% Case for no treatment
[Tx_uncon_CD,X_uncon_CD] = ode15s(@(t,x) stateEqCD(t,x,uCD,Tu), Tu, ...
    initx(1:4), options);
[Tx_uncon_HFD,X_uncon_HFD] = ode15s(@(t,x) stateEqHFD(t,x,uHFD,Tu), Tu, ...
    initx(5:end), options);

Tumor=ceil(X_uncon_CD(:,1)+ X_uncon_CD(:,2));
index_CD=(find(Tumor>tumor_threshold))-1;
t0_treatment_CD=(Tx_uncon_CD(min(index_CD)));
Tumor(min(index_CD));

% Time discretization for CD
t0_treatment=t0_treatment_CD; % treatment starts at t=t0_treatment
Tu_treatment=t0_treatment:0.1:tf;

Tumor=ceil(X_uncon_HFD(:,1)+ X_uncon_HFD(:,2));
index_HFD=find(Tumor>tumor_threshold)-1;
t0_treatment_HFD= (Tx_uncon_HFD(min(index_HFD)));
Tumor(min(index_HFD));

%Preperation for treatment
[~,X_uncon_treat_CD] = ode15s(@(t,x) stateEqCD(t,x,uCD,Tu), [t0:.0001:t0_treatment_CD], ...
    initx(1:4), options);
[~,X_uncon_treat_HFD] = ode15s(@(t,x) stateEqHFD(t,x,uHFD,Tu), [t0:.0001:t0_treatment_HFD], ...
    initx(5:end), options);

% Initial condition for treatment
initx_treatment_CD = X_uncon_treat_CD(end,:);
initx_treatment_HFD = X_uncon_treat_HFD(end,:);

%% Adaptive treatment

% Set the treatment for CD
uCD_adapt=0*ones(size(Tu_treatment))';
for ii=1:2:6%6 % 10
    uCD_adapt((ii-1)*21+1:((ii*21)))=M2*1; % 21 11
end

% Solve the DE
[Tx,X] = ode15s(@(t,x) stateEqCD(t,x,uCD_adapt,Tu_treatment), Tu_treatment, ...
    initx_treatment_CD, options);

XCD=X; TxCD=Tx;

% Time discretization for HFD
t0_treatment=t0_treatment_HFD; % treatment starts at t=t0_treatment
Tu_treatment=t0_treatment:0.1:tf;

% Set the treatment for HFD
uHFD_adapt=0*ones(size(Tu_treatment))';
for ii=1:2:6%6 %6 11
    uHFD_adapt((ii-1)*21+1:((ii*21)))=M2*1; % 21 11
end

% Solve the DE
[Tx,X] = ode15s(@(t,x) stateEqHFD(t,x,uHFD_adapt,Tu_treatment), Tu_treatment, ...
    initx_treatment_HFD, options);

XHFD=X; TxHFD=Tx;

%%

figure(57)
colororder({'b','k'})
yyaxis left
f1=plot(Tx_uncon_CD,X_uncon_CD(:,1)+X_uncon_CD(:,2),'-','LineWidth',4);
hold on
f2=plot(TxCD,XCD(:,1)+XCD(:,2),':','LineWidth',4);
grid on
xlabel('t')
ylabel('S(t)+R(t)')
xlim([0,tf])
ylim([0,2100])
yyaxis right
f3=stairs(TxCD,uCD_adapt,'--','LineWidth',2);
ylabel('u(t)')
legend([f1,f2,f3],{'u=0','u \neq 0','Control func.'},'Location','NorthWest')
title('CD');
xline(t0_treatment_CD,'--k','HandleVisibility','off')

figure(58)
colororder({'b','k'})
yyaxis left
f1=plot(Tx_uncon_HFD,X_uncon_HFD(:,1)+X_uncon_HFD(:,2),'-','LineWidth',4);
hold on
f2=plot(TxHFD,XHFD(:,1)+XHFD(:,2),':','LineWidth',4);
grid on
legend('u_{HFD}=0','u_{HFD} \neq 0','Location','SouthEast')
xlabel('t')
ylabel('S(t)+R(t)')
xlim([0,tf])
ylim([0,2100])
yyaxis right
f3=stairs(TxHFD,uHFD_adapt,'-.','LineWidth',2);
ylabel('u(t)')
legend([f1,f2,f3],{'u=0','u \neq 0','Control func.'},'Location','NorthWest')
title('HFD');
xline(t0_treatment_HFD,'-.k','HandleVisibility','off')

end
