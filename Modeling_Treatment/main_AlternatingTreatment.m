%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to simulate alternating treatment
%
% Author: Tuğba Akman Date: Jan 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main_AlternatingTreatment()

close all
clear
clc
format long

set(0, 'defaultaxesfontsize',16)

global Eq


%Initial conditions
T0=1;
S0=1;
R0=T0-S0;
E0_CD=1040.35/5.94;
E0_HFD=7685.87/5.94;
F0_CD = 49.923;
F0_HFD = 368.82;

% Parameters
Eq.mu = 5.94;
Eq.m1 = 5e-4;
Eq.m2CD = 1/(F0_HFD);
Eq.m2HFD = 1/(F0_HFD);
Eq.n = 1;
Eq.k1 = 0.586967;
Eq.a1 = 59.0927;
Eq.r = 20.8391;
Eq.alpha = 2.21427e-05;
Eq.pCD=1;
Eq.pHFD=1;
Eq.K=1/Eq.m1;
Eq.m=1;
Eq.betaCD=1;
Eq.betaHFD=1;
Eq.l=10;
Eq.Sbound=0;
Eq.k2CD = 0.045;
Eq.k2HFD= 0.045;
Eq.k3CD=0.25*Eq.k1;
Eq.k3HFD=0.25*Eq.k1;
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
t0_treatment_CD=(Tx_uncon_CD(min(index_CD)))
Tumor(min(index_CD))

% Time discretization for CD
t0_treatment=t0_treatment_CD; % treatment starts at t=t0_treatment
Tu_treatment=t0_treatment:0.1:tf;

Tumor=ceil(X_uncon_HFD(:,1)+ X_uncon_HFD(:,2));
index_HFD=find(Tumor>tumor_threshold)-1;
t0_treatment_HFD= (Tx_uncon_HFD(min(index_HFD)))
Tumor(min(index_HFD))

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
for ii=1:2:5
    uCD_adapt((ii-1)*21+1:((ii*21)))=M2*1;
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
for ii=1:2:6
    uHFD_adapt((ii-1)*22+1:((ii*22)))=M2*1;
end

% Solve the DE
[Tx,X] = ode15s(@(t,x) stateEqHFD(t,x,uHFD_adapt,Tu_treatment), Tu_treatment, ...
    initx_treatment_HFD, options);

XHFD=X; TxHFD=Tx;

%%
figure(51)
plot(TxHFD,M2*ones(size(TxHFD)),'-r','LineWidth',4)
hold on
stairs(TxCD,uCD_adapt,'c-','LineWidth',4)
hold on
stairs(TxHFD,uHFD_adapt,'b:.','LineWidth',4)
title('Control function')
xlabel('t')
ylabel('u(t)')
legend('Upper bound','CD','HFD','Location','SouthEast')
xlim([10,tf])
grid on


figure(52)
subplot(2,4,1)
plot(Tx_uncon_CD,X_uncon_CD(:,1),'r-','LineWidth',4)
hold on
plot(TxCD,XCD(:,1),'b:.','LineWidth',4)
grid on
xlabel('t')
ylabel('S(t)')
title('CD')
ylim([0,2100])
xlim([0,tf])

subplot(2,4,2)
plot(Tx_uncon_CD,X_uncon_CD(:,2),'r-','LineWidth',4)
hold on
plot(TxCD,XCD(:,2),'b:.','LineWidth',4)
grid on
xlabel('t')
ylabel('R(t)')
title('CD')
ylim([0,2100])
xlim([0,tf])

subplot(2,4,3)
plot(Tx_uncon_CD,X_uncon_CD(:,3),'r-','LineWidth',4)
hold on
plot(TxCD,XCD(:,3),'b:.','LineWidth',4)
grid on
legend('u_{CD}=0','u_{CD} \neq 0','Location','NorthEast')
xlabel('t')
ylabel('E(t)')
title('CD')
ylim([0,1500])
xlim([0,tf])

subplot(2,4,4)
plot(Tx_uncon_CD,X_uncon_CD(:,4),'r-','LineWidth',4)
hold on
plot(TxCD,XCD(:,4),'b:.','LineWidth',4)
grid on
xlabel('t')
ylabel('F(t)')
title('CD')
ylim([0,400])
xlim([0,tf])

subplot(2,4,5)
plot(Tx_uncon_HFD,X_uncon_HFD(:,1),'r-','LineWidth',4)
hold on
plot(TxHFD,XHFD(:,1),'b:.','LineWidth',4)
grid on
xlabel('t')
ylabel('S(t)')
title('HFD')
ylim([0,2100])
xlim([0,tf])

subplot(2,4,6)
plot(Tx_uncon_HFD,X_uncon_HFD(:,2),'r-','LineWidth',4)
hold on
plot(TxHFD,XHFD(:,2),'b:.','LineWidth',4)
grid on
xlabel('t')
ylabel('R(t)')
title('HFD')
ylim([0,2100])
xlim([0,tf])

subplot(2,4,7)
plot(Tx_uncon_HFD,X_uncon_HFD(:,3),'r-','LineWidth',4)
hold on
plot(TxHFD,XHFD(:,3),'b:.','LineWidth',4)
grid on
legend('u_{HFD}=0','u_{HFD} \neq 0','Location','NorthEast')
xlabel('t')
ylabel('E(t)')
title('HFD')
ylim([0,1500])
xlim([0,tf])

subplot(2,4,8)
plot(Tx_uncon_HFD,X_uncon_HFD(:,4),'r-','LineWidth',4)
hold on
plot(TxHFD,XHFD(:,4),'b:.','LineWidth',4)
grid on
xlabel('t')
ylabel('F(t)')
title('HFD')
ylim([0,400])
xlim([0,tf])

figure(55)
plot(Tx_uncon_CD,X_uncon_CD(:,1)+X_uncon_CD(:,2),'r-','LineWidth',4)
hold on
plot(TxCD,XCD(:,1)+XCD(:,2),'b:.','LineWidth',4)
grid on
legend('u_{CD}=0','u_{CD} \neq 0','Location','SouthEast')
xlabel('t')
ylabel('S(t)+R(t)')
xlim([0,tf])
ylim([0,2100])
title('CD');

figure(56)
plot(Tx_uncon_HFD,X_uncon_HFD(:,1)+X_uncon_HFD(:,2),'r-','LineWidth',4)
hold on
plot(TxHFD,XHFD(:,1)+XHFD(:,2),'b:.','LineWidth',4)
grid on
legend('u_{HFD}=0','u_{HFD} \neq 0','Location','SouthEast')
xlabel('t')
ylabel('S(t)+R(t)')
xlim([0,tf])
ylim([0,2100])
title('HFD');

figure(57)
yyaxis left
f1=plot(Tx_uncon_CD,X_uncon_CD(:,1)+X_uncon_CD(:,2),'-','LineWidth',4);
hold on
f2=plot(TxCD,XCD(:,1)+XCD(:,2),':.','LineWidth',4);
grid on
xlabel('t')
ylabel('S(t)+R(t)')
xlim([0,tf])
ylim([0,2100])
yyaxis right
f3=stairs(TxCD,uCD_adapt,'-','LineWidth',2);
ylabel('u(t)')
legend([f1,f2,f3],{'u=0','u \neq 0','Control func.'},'Location','NorthWest')
title('CD');

figure(58)
yyaxis left
f1=plot(Tx_uncon_HFD,X_uncon_HFD(:,1)+X_uncon_HFD(:,2),'-','LineWidth',4);
hold on
f2=plot(TxHFD,XHFD(:,1)+XHFD(:,2),':.','LineWidth',4);
grid on
legend('u_{HFD}=0','u_{HFD} \neq 0','Location','SouthEast')
xlabel('t')
ylabel('S(t)+R(t)')
xlim([0,tf])
ylim([0,2100])
yyaxis right
f3=stairs(TxHFD,uHFD_adapt,'-','LineWidth',2);
ylabel('u(t)')
legend([f1,f2,f3],{'u=0','u \neq 0','Control func.'},'Location','NorthWest')
title('HFD');

end
