%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to simulate optimal anti-hormonal treatment
%
% Author: Tuğba Akman Date: July 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data are taken from
% Obesity-Activated Adipose-Derived Stromal Cells Promote
% Breast Cancer Growth and Invasion
% Neoplasia (2018) 20, 1161-1174
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main_OCP()

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
Eq.a2=20;
Eq.a3=1;

% Weights in the cost functional
Eq.weight1 = 1; % omega_S for CD
Eq.weight2 = 1; % omega_R for CD
Eq.weight3 = 1; % omega_U for CD
Eq.weight4 = 1; % omega_S for HFD 
Eq.weight5 = 1; % omega_R for HFD
Eq.weight6 = 1; % omega_U for HFD

% Control constraints
M1=0;    % Lower bound for control
M2=0.99; % Upper bound for control

% Time discretization
t0=0;                   % initial time
tf = 25;                % final time
Interval=100;           % Number of subintervals
dt = (tf-t0)/Interval; % Temporal increment
Tu=t0:(1*dt):tf;        % ınterpolation points

% tumor_threshold to start treatment
tumor_threshold = ceil(Eq.K*0.25);

initx=[S0 R0 E0_CD F0_CD S0 R0 E0_HFD F0_HFD]'; % IC for state
initlambda=[0,0,0,0,0,0,0,0]';                  % FC for adjoint

%Allocate for uncontrolled case - MUST BE EQUAL TO ZERO!!!!
uCD=zeros(size(Tu))';
uHFD=zeros(size(Tu))';

options = odeset('RelTol',1e-6,'AbsTol',1e-6);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Uncontrolled case - No treatment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Tx_uncon_CD,X_uncon_CD] = ode15s(@(t,x) stateEqCD(t,x,uCD,Tu), Tu, ...
    initx(1:4), options);
[Tx_uncon_HFD,X_uncon_HFD] = ode15s(@(t,x) stateEqHFD(t,x,uHFD,Tu), Tu, ...
    initx(5:end), options);

figure;
plot(Tx_uncon_CD,X_uncon_CD(:,1))
%pause
% Time to start treatment for CD
Tumor=ceil(X_uncon_CD(:,1)+ X_uncon_CD(:,2))
index_CD=(find(Tumor>tumor_threshold))-1;
t0_treatment_CD=(Tx_uncon_CD(min(index_CD)))
Tumor(min(index_CD));

% Time discretization for treatment - CD
t0_treatment=t0_treatment_CD;           % treatment starts at t=t0_treatment
Interval_treatment=100;                 % Number of subintervals
dt_treatment = ( tf-t0_treatment)/Interval_treatment; % Temporal increment
Tu_treatment=t0_treatment:(1*dt_treatment):tf;        % interpolation points

% Time to start treatment for HFD
Tumor=ceil(X_uncon_HFD(:,1)+ X_uncon_HFD(:,2));
index_HFD=find(Tumor>tumor_threshold)-1;
t0_treatment_HFD= (Tx_uncon_HFD(min(index_HFD)))
Tumor(min(index_HFD));
pause
%Preperation for treatment
[Tx_uncon_treat_CD,X_uncon_treat_CD] = ode15s(@(t,x) stateEqCD(t,x,uCD,Tu), [t0:.0001:t0_treatment_CD], ...
    initx(1:4), options);
[Tx_uncon_treat_HFD,X_uncon_treat_HFD] = ode15s(@(t,x) stateEqHFD(t,x,uHFD,Tu), [t0:.0001:t0_treatment_HFD], ...
    initx(5:end), options);

% Initial condition for treatment
initx_treatment_CD = X_uncon_treat_CD(end,:);
initx_treatment_HFD = X_uncon_treat_HFD(end,:);

% Allocate for treatment
uCD=0*ones(size(Tu_treatment))'+0.95;
uHFD=0*ones(size(Tu_treatment))'+0.95;

theta_convex=linspace(0.01,0.99,100);

lambda1=zeros(size(Tu_treatment));
lambda2=zeros(size(Tu_treatment));
lambda3=zeros(size(Tu_treatment));
lambda4=zeros(size(Tu_treatment));
x1=zeros(size(Tu_treatment));
x2=zeros(size(Tu_treatment));
x3=zeros(size(Tu_treatment));
x4=zeros(size(Tu_treatment));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solve OCP for CD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nn=1; test = 11;
while(test > 1e-5)

    olduCD = uCD;
    oldx1=x1;
    oldx2=x2;
    oldx3=x3;
    oldx4=x4;
    oldlambda1=lambda1;
    oldlambda2=lambda2;
    oldlambda3=lambda3;
    oldlambda4=lambda4;
    
    % solve the state eqn forward
    [Tx,X] = ode15s(@(t,x) stateEqCD(t,x,uCD,Tu_treatment), Tu_treatment, ...
        initx_treatment_CD, options);

    % solve the adjoint eqn backward
    x1 = X(:,1); x2 = X(:,2); x3 = X(:,3); x4 = X(:,4);
    [Tlambda,Lambda] = ode15s(@(t,lambda) costateEqCD(t,lambda,uCD,Tu_treatment,x1,x2,x3,x4,Tx), ...
        flip(Tu_treatment), initlambda(1:4), options);

    lambda1 = Lambda(:,1); lambda2 = Lambda(:,2); lambda3 = Lambda(:,3);
    lambda4 = Lambda(:,4);

    % Interpolation
    lambda1 = interp1(Tlambda,lambda1,Tx);
    lambda2 = interp1(Tlambda,lambda2,Tx);
    lambda3 = interp1(Tlambda,lambda3,Tx);
    lambda4 = interp1(Tlambda,lambda4,Tx);

    % Optimality condition
    optCond1 = lambda3.*x4.*Eq.pCD.*Eq.r/Eq.weight3;

    % Project the control
    u1CD = min(M2, max(M1, optCond1));

    % Update the control
    for dd=1:1:length(theta_convex)

        thetaa=theta_convex(dd);
        uCD_convex_for_inside(:,dd) = (thetaa)*u1CD + (1-thetaa)*olduCD;       
        JCD_for_inside(dd)=cost(Tx,X,uCD_convex_for_inside(:,dd));
        
    end

    [JCD_min_for_inside_val, pos_min_for_inside] = min(JCD_for_inside);

    uCD_convex_for_inside(:,pos_min_for_inside);

    uCD = uCD_convex_for_inside(:,pos_min_for_inside);

    %Use same interpolation  points
    Tu_treatment=Tx;

    % Stopping criteria on relative error
    temp1 = norm((olduCD - uCD))/norm((uCD));
    temp2 = norm((oldx1 - x1))/norm(x1);
    temp3 = norm((oldx2 - x2))/norm(x2);
    temp4 = norm((oldx3 - x3))/norm(x3);
    temp5 = norm((oldx4 - x4))/norm(x4);

    temp6 = norm((oldlambda1 - lambda1))/norm((lambda1));
    temp7 = norm((oldlambda2 - lambda2))/norm((lambda2));
    temp8 = norm((oldlambda3 - lambda3))/norm((lambda3));
    temp9 = norm((oldlambda4 - lambda4))/norm((lambda4));

    test = max(temp1, max(temp2, max( temp3, max( temp4, max(temp5,...
        max(temp6, max(temp7, max(temp8, temp9))))))))

    %test=temp1;

    XCD=X; TxCD=Tx;

    JCD_new(nn)=cost(TxCD,XCD,uCD); 

    nn=nn+1;

end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solve OCP for HFD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time discretization for treatment - HFD
t0_treatment=t0_treatment_HFD; % treatment starts at t=t0_treatment
Interval_treatment=100; % Number of subintervals
%Nf_treatment=Interval_treatment+1; %  index of the last time step
dt_treatment = ( tf-t0_treatment)/Interval_treatment; % Temporal increment
%t_treatment = t0_treatment:dt_treatment:tf; % Discrete time
Tu_treatment=t0_treatment:(1*dt_treatment):tf; % ınterpolation points


lambda1=zeros(size(Tu_treatment));
lambda2=zeros(size(Tu_treatment));
lambda3=zeros(size(Tu_treatment));
lambda4=zeros(size(Tu_treatment));
x1=zeros(size(Tu_treatment));
x2=zeros(size(Tu_treatment));
x3=zeros(size(Tu_treatment));
x4=zeros(size(Tu_treatment));

test = 1; nn=1;

while(test > 1e-5)

    olduHFD = uHFD;
    oldx1=x1;
    oldx2=x2;
    oldx3=x3;
    oldx4=x4;
    oldlambda1=lambda1;
    oldlambda2=lambda2;
    oldlambda3=lambda3;
    oldlambda4=lambda4;   

    % solve the state eqn forward
    [Tx,X] = ode15s(@(t,x) stateEqHFD(t,x,uHFD,Tu_treatment), Tu_treatment, ...
        initx_treatment_HFD, options);

    % solve the adjoint eqn backward
    x1 = X(:,1); x2 = X(:,2); x3 = X(:,3); x4 = X(:,4);
    [Tlambda,Lambda] = ode15s(@(t,lambda) costateEqHFD(t,lambda,uHFD,Tu_treatment,x1,x2,x3,x4,Tx), ...
        flip(Tu_treatment), initlambda(5:end), options);

    lambda1 = Lambda(:,1); lambda2 = Lambda(:,2); lambda3 = Lambda(:,3);
    lambda4 = Lambda(:,4);

    % Interpolation
    lambda1 = interp1(Tlambda,lambda1,Tx);
    lambda2 = interp1(Tlambda,lambda2,Tx);
    lambda3 = interp1(Tlambda,lambda3,Tx);
    lambda4 = interp1(Tlambda,lambda4,Tx);
   
    % Optimality condition
    optCond1 = lambda3.*x4.*Eq.pHFD.*Eq.r/Eq.weight6;

    % Project the control
    u1HFD = min(M2, max(M1, optCond1));
    
    % Update the control
    for dd=1:1:length(theta_convex)

        thetaa=theta_convex(dd);
        uHFD_convex_for_inside(:,dd) = (thetaa)*u1HFD + (1-thetaa)*olduHFD;
        JHFD_for_inside(dd)=cost(Tx,X,uHFD_convex_for_inside(:,dd));
        
    end

    [JHFD_min_for_inside_val, pos_min_for_inside] = min(JHFD_for_inside);

    uHFD_convex_for_inside(:,pos_min_for_inside);

    uHFD = uHFD_convex_for_inside(:,pos_min_for_inside);
    
    %Use same interpolation  points
    Tu_treatment=Tx;

    % Stopping criteria on relative error
    temp1 = norm((olduHFD - uHFD))/norm((uHFD));
    temp2 = norm((oldx1 - x1))/norm(x1);
    temp3 = norm((oldx2 - x2))/norm(x2);
    temp4 = norm((oldx3 - x3))/norm(x3);
    temp5 = norm((oldx4 - x4))/norm(x4);

    temp6 = norm((oldlambda1 - lambda1))/norm((lambda1));
    temp7 = norm((oldlambda2 - lambda2))/norm((lambda2));
    temp8 = norm((oldlambda3 - lambda3))/norm((lambda3));
    temp9 = norm((oldlambda4 - lambda4))/norm((lambda4));

    test = max(temp1, max(temp2, max( temp3, max( temp4, max(temp5,...
        max(temp6, max(temp7, max(temp8, temp9))))))))
    
    XHFD=X; TxHFD=Tx;

    JHFD_new(nn)=cost(TxHFD,XHFD,uHFD);

    nn=nn+1;

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
plot(TxHFD,M2*ones(size(TxHFD)),'-m','LineWidth',4)
hold on
plot(TxCD,uCD,'k--','LineWidth',2)
hold on
plot(TxHFD,uHFD,'k-.','LineWidth',2)
title('Control function')
xlabel('t')
ylabel('u(t)')
legend('Maximum treatment','CD','HFD','Location','SouthWest')
xlim([10,tf])
xline(t0_treatment_CD,'--k','HandleVisibility','off')
hold on
xline(t0_treatment_HFD,'-.k','HandleVisibility','off')
grid on

XCD=[X_uncon_treat_CD; XCD];
TxCD=[Tx_uncon_treat_CD;TxCD];
XHFD=[X_uncon_treat_HFD; XHFD];
TxHFD=[Tx_uncon_treat_HFD;TxHFD];

figure(330)
plot(Tx_uncon_CD,X_uncon_CD(:,1)+X_uncon_CD(:,2),'b-','LineWidth',4)
hold on
plot(TxCD,XCD(:,1)+XCD(:,2),'b:','LineWidth',4)
grid on
legend('u_{CD}=0','u_{CD} \neq 0','Location','NorthWest')
xlabel('t')
ylabel('S(t)+R(t)')
xlim([0,tf])
ylim([0,2100])
title('CD');
xline(t0_treatment_CD,'--k','HandleVisibility','off')

figure(331)
plot(Tx_uncon_HFD,X_uncon_HFD(:,1)+X_uncon_HFD(:,2),'b-','LineWidth',4)
hold on
plot(TxHFD,XHFD(:,1)+XHFD(:,2),'b:','LineWidth',4)
grid on
legend('u_{HFD}=0','u_{HFD} \neq 0','Location','NorthWest')
xlabel('t')
ylabel('S(t)+R(t)')
xlim([0,tf])
ylim([0,2100])
title('HFD');
xline(t0_treatment_HFD,'--k','HandleVisibility','off')

figure(21)
subplot(2,4,1)
plot(Tx_uncon_CD,X_uncon_CD(:,1),'b-','LineWidth',4)
hold on
plot(TxCD,XCD(:,1),'b:','LineWidth',4)
grid on
xlabel('t')
ylabel('S(t)')
title('CD')
ylim([0,2100])
xlim([0,tf])
xline(t0_treatment_CD,'--k','HandleVisibility','off')

subplot(2,4,2)
plot(Tx_uncon_CD,X_uncon_CD(:,2),'b-','LineWidth',4)
hold on
plot(TxCD,XCD(:,2),'b:','LineWidth',4)
grid on
xlabel('t')
ylabel('R(t)')
title('CD')
ylim([0,2100])
xlim([0,tf])
legend('u_{CD}=0','u_{CD} \neq 0','Location','NorthWest')
xline(t0_treatment_CD,'--k','HandleVisibility','off')

subplot(2,4,3)
plot(Tx_uncon_CD,X_uncon_CD(:,3),'b-','LineWidth',4)
hold on
plot(TxCD,XCD(:,3),'b:','LineWidth',4)
grid on
xlabel('t')
ylabel('E(t)')
title('CD')
ylim([0,1500])
xlim([0,tf])
xline(t0_treatment_CD,'--k','HandleVisibility','off')

subplot(2,4,4)
plot(Tx_uncon_CD,X_uncon_CD(:,4),'b-','LineWidth',4)
hold on
plot(TxCD,XCD(:,4),'b:','LineWidth',4)
grid on
xlabel('t')
ylabel('F(t)')
title('CD')
ylim([0,1100])
xlim([0,tf])
xline(t0_treatment_CD,'--k','HandleVisibility','off')

subplot(2,4,5)
plot(Tx_uncon_HFD,X_uncon_HFD(:,1),'b-','LineWidth',4)
hold on
plot(TxHFD,XHFD(:,1),'b:','LineWidth',4)
grid on
xlabel('t')
ylabel('S(t)')
title('HFD')
ylim([0,2100])
xlim([0,tf])
xline(t0_treatment_HFD,'--k','HandleVisibility','off')

subplot(2,4,6)
plot(Tx_uncon_HFD,X_uncon_HFD(:,2),'b-','LineWidth',4)
hold on
plot(TxHFD,XHFD(:,2),'b:','LineWidth',4)
grid on
xlabel('t')
ylabel('R(t)')
title('HFD')
ylim([0,2100])
xlim([0,tf])
legend('u_{HFD}=0','u_{HFD} \neq 0','Location','NorthWest')
xline(t0_treatment_HFD,'--k','HandleVisibility','off')

subplot(2,4,7)
plot(Tx_uncon_HFD,X_uncon_HFD(:,3),'b-','LineWidth',4)
hold on
plot(TxHFD,XHFD(:,3),'b:','LineWidth',4)
grid on
xlabel('t')
ylabel('E(t)')
title('HFD')
ylim([0,1500])
xlim([0,tf])
xline(t0_treatment_HFD,'--k','HandleVisibility','off')

subplot(2,4,8)
plot(Tx_uncon_HFD,X_uncon_HFD(:,4),'b-','LineWidth',4)
hold on
plot(TxHFD,XHFD(:,4),'b:','LineWidth',4)
grid on
xlabel('t')
ylabel('F(t)')
title('HFD')
ylim([0,1100])
xlim([0,tf])
xline(t0_treatment_HFD,'--k','HandleVisibility','off')
end

%Sub-function

function J=cost(Tx,X,u)
global Eq

J=(Eq.weight1*trapz(Tx,X(:,1)-Eq.Sbound) + Eq.weight2*trapz(Tx,X(:,2)) + Eq.weight3*0.5*trapz(Tx,u.^2));
end

