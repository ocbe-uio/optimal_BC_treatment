%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to simulate extended model
%
% Author: Tugba Akman Date: Jan 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data are taken from
% Obesity-Activated Adipose-Derived Stromal Cells Promote
% Breast Cancer Growth and Invasion
% Neoplasia (2018) 20, 1161–1174
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main_ExtendedModel

clear
close all
clc
format long

set(0, 'defaultaxesfontsize',16)

global tumor_CD_tdata tumor_CD_ydata tumor_HFD_tdata tumor_HFD_ydata...
    NumbOfASC_CD_ydata NumbOfASC_HFD_ydata...
    diam_CD diam_HFD Fat_CD_ydata Fat_CD_tdata Fat_HFD_ydata Fat_HFD_tdata...
    E0_CD E0_HFD T0 F0_CD F0_HFD mu K k3CD k3HFD m betaCD betaHFD m2CD ...
    m2HFD a1 n alpha k1 k2CD k2HFD dCD dHFD pCD pHFD r m1 S0 R0 a2 l a3 cc

%Data
tumor_CD_tdata = [10, 13, 15]';
tumor_CD_ydata = [138.92, 324.7, 480.76]';
tumor_HFD_tdata = [10, 13, 15]';
tumor_HFD_ydata = [326.8, 702.025, 1114.833]';
Fat_CD_tdata = 15;
Fat_HFD_tdata = 15;

%Formula for # of ASCs
diam_CD = 0.1;
diam_HFD = 0.2;
radius_CD=diam_CD/2;
radius_HFD=diam_HFD/2;
NumbOfASC_CD_ydata = ((tumor_CD_ydata(3))^(2/3))*(143.18/diam_CD); 
NumbOfASC_HFD_ydata = ((tumor_HFD_ydata(3))^(2/3))*(159.09/diam_HFD); 

% Amount of fat
Fat_CD_ydata = NumbOfASC_CD_ydata*(4/3)*pi*(radius_CD^3); 
Fat_HFD_ydata = NumbOfASC_HFD_ydata*(4/3)*pi*(radius_HFD^3); 

%Initial conditions
T0=1;
S0=1;
R0=T0-S0;
E0_CD=1040.35/5.94;
E0_HFD=7685.87/5.94;
F0_CD = 49.923;
F0_HFD = 368.82;

% Time discretization
dt = 0.01;
t0CD=13.75;
t0HFD=11.5;
tf=25;
tforwardCD = 0:dt:t0CD; % t mesh for CD
tforwardHFD = 0:dt:t0HFD; % t mesh for HFD
tforwardCD_treat = t0CD:dt:tf; % t mesh CD during treatment
tforwardHFD_treat = t0HFD:dt:tf; % t mesh HFD during treatment

% Parameters
mu = 5.94;
m1 = 5e-4;
m2CD = 1/(F0_HFD); 
m2HFD = 1/(F0_HFD);
n = 1;
k1 = 0.586967;
a1 = 59.0927;
r = 20.8391;
alpha = 2.21427e-05;
cc=1;
l=10;
dCD=0;
dHFD=0;
K=1/m1;
m=1;
betaCD=1;
betaHFD=1;
k2CD = 0.045;
k2HFD= 0.045;
a2=10;
a3=10;
k3CD=0.25*k1;
k3HFD=0.25*k1;

pCD_test=0.5*[2,1/20,1/40,1/50,1/500];
pHFD_test=0.5*[2,1/20,1/40,1/50,1/500];

params=[];

options = odeset('RelTol',1e-6,'AbsTol',1e-6);

% IC for phase2
initial_cond = [S0 R0 E0_CD F0_CD S0 R0 E0_HFD F0_HFD];

for i=1:length(pCD_test)

    pCD=1;
    pHFD=1;
    [~,X_uncon_treat_CD] = ode15s(@(t,y)Model_CD(y,params),tforwardCD,initial_cond(1:4),options);
    [~,X_uncon_treat_HFD] = ode15s(@(t,y)Model_HFD(y,params),tforwardHFD,initial_cond(5:end),options);

    pCD = pCD_test(i);
    pHFD = pHFD_test(i);

    %Solve the model after data-fitting
    [~, y_r_phase2_CD] = ode15s(@(t,y)Model_CD(y,params),tforwardCD_treat,X_uncon_treat_CD(end,:),options);
    [~, y_r_phase2_HFD] = ode15s(@(t,y)Model_HFD(y,params),tforwardHFD_treat,X_uncon_treat_HFD(end,:),options);

    out_big_CD(:,:,i) =  y_r_phase2_CD;
    out_big_HFD(:,:,i) =  y_r_phase2_HFD;
end


pCD=1;
pHFD=1;
[Tx_uncon_treat_CD,X_uncon_treat_CD] = ode15s(@(t,y)Model_CD(y,params),tforwardCD,initial_cond(1:4),options);
[Tx_uncon_treat_HFD,X_uncon_treat_HFD] = ode15s(@(t,y)Model_HFD(y,params),tforwardHFD,initial_cond(5:end),options);



%%
colormap lines

figure(113)
colors = jet(5);
%subplot(1,2,1)
plot(tforwardCD_treat,out_big_CD(:,1,1) + m*out_big_CD(:,2,1),'r','LineWidth',2)
hold on
plot(tforwardCD_treat,out_big_CD(:,1,2) + m*out_big_CD(:,2,2),'c-','LineWidth',4)
hold on
plot(tforwardCD_treat,out_big_CD(:,1,3)+m*out_big_CD(:,2,3),'-.','color',colors(3,:),'LineWidth',5);
hold on
plot(tforwardCD_treat,out_big_CD(:,1,4)+m*out_big_CD(:,2,4),'m--','LineWidth',4);
hold on
plot(tforwardCD_treat,out_big_CD(:,1,5)+m*out_big_CD(:,2,5),':','color',colors(1,:),'LineWidth',4);
hold on
plot(tforwardCD_treat,out_big_CD(:,1,5)+m*out_big_CD(:,2,5),':','color',colors(1,:),'LineWidth',4);
hold on
plot(Tx_uncon_treat_CD,X_uncon_treat_CD(:,1),'r','LineWidth',2);
xlabel('t')
ylim([0,2100])
title('CD');
ylabel('S(t)+R(t)')
legend(strcat('p= ',num2str(pCD_test')),'Location','Best','Orientation','vertical')
xlim([0,tf])
grid on

colors = jet(5);
figure(114)
plot(tforwardHFD_treat,out_big_HFD(:,1,1) + m*out_big_HFD(:,2,1),'r','LineWidth',2)
hold on
plot(tforwardHFD_treat,out_big_HFD(:,1,2) + m*out_big_HFD(:,2,2),'c-','LineWidth',4)
hold on
plot(tforwardHFD_treat,out_big_HFD(:,1,3)+m*out_big_HFD(:,2,3),'-.','color',colors(3,:),'LineWidth',5);
hold on
plot(tforwardHFD_treat,out_big_HFD(:,1,4)+m*out_big_HFD(:,2,4),'m--','LineWidth',4);
hold on
plot(tforwardHFD_treat,out_big_HFD(:,1,5)+m*out_big_HFD(:,2,5),':','color',colors(1,:),'LineWidth',4);
hold on
plot(Tx_uncon_treat_HFD,X_uncon_treat_HFD(:,1),'r','LineWidth',2);
xlabel('t')
ylabel('S(t)+R(t)')
ylim([0,2100])
title('HFD');
xlim([0,tf])
grid on


figure(111)
colors = jet(5);
subplot(2,4,1)
plot(tforwardCD_treat,out_big_CD(:,1,1),'color',colors(5,:),'LineWidth',2)
hold on
plot(tforwardCD_treat,out_big_CD(:,1,2),'c-','LineWidth',4)
hold on
plot(tforwardCD_treat,out_big_CD(:,1,3),'-.','color',colors(3,:),'LineWidth',5);
hold on
plot(tforwardCD_treat,out_big_CD(:,1,4),'m--','LineWidth',4);
hold on
plot(tforwardCD_treat,out_big_CD(:,1,5),':','color',colors(1,:),'LineWidth',4);
hold on
plot(Tx_uncon_treat_CD,X_uncon_treat_CD(:,1),'color',colors(5,:),'LineWidth',2);
xlabel('t')
ylabel('S(t)')
title('CD');
ylim([0,2100])
xlim([0,tf])
grid on

subplot(2,4,2)
colors = jet(5);
plot(tforwardCD_treat,out_big_CD(:,2,1),'color',colors(5,:),'LineWidth',2)
hold on
plot(tforwardCD_treat,out_big_CD(:,2,2),'c-','LineWidth',4)
hold on
plot(tforwardCD_treat,out_big_CD(:,2,3),'-.','color',colors(3,:),'LineWidth',5);
hold on
plot(tforwardCD_treat,out_big_CD(:,2,4),'m--','LineWidth',4);
hold on
plot(tforwardCD_treat,out_big_CD(:,2,5),':','color',colors(1,:),'LineWidth',4);
hold on
plot(Tx_uncon_treat_CD,X_uncon_treat_CD(:,2),'color',colors(5,:),'LineWidth',2);
xlabel('t')
ylim([0,2100])
ylabel('R(t)')
title('CD');
xlim([0,tf])
grid on

subplot(2,4,3)
colors = jet(5);
plot(tforwardCD_treat,out_big_CD(:,3,1),'color',colors(5,:),'LineWidth',2)
hold on
plot(tforwardCD_treat,out_big_CD(:,3,2),'c-','LineWidth',4)
hold on
plot(tforwardCD_treat,out_big_CD(:,3,3),'-.','color',colors(3,:),'LineWidth',5);
hold on
plot(tforwardCD_treat,out_big_CD(:,3,4),'m--','LineWidth',4);
hold on
plot(tforwardCD_treat,out_big_CD(:,3,5),':','color',colors(1,:),'LineWidth',4);
hold on
plot(Tx_uncon_treat_CD,X_uncon_treat_CD(:,3),'color',colors(5,:),'LineWidth',2);
xlabel('t')
ylim([0,1500])
ylabel('E(t)')
title('CD');
xlim([0,tf])
grid on
legend(strcat('p_{CD}= ',num2str(pCD_test')),'Location','NorthEast','Orientation','vertical')

subplot(2,4,4)
colors = jet(5);
plot(tforwardCD_treat,out_big_CD(:,4,1),'color',colors(5,:),'LineWidth',2)
hold on
plot(tforwardCD_treat,out_big_CD(:,4,2),'c-','LineWidth',4)
hold on
plot(tforwardCD_treat,out_big_CD(:,4,3),'-.','color',colors(3,:),'LineWidth',5);
hold on
plot(tforwardCD_treat,out_big_CD(:,4,4),'m--','LineWidth',4);
hold on
plot(tforwardCD_treat,out_big_CD(:,4,5),':','color',colors(1,:),'LineWidth',4);
hold on
plot(Tx_uncon_treat_CD,X_uncon_treat_CD(:,4),'color',colors(5,:),'LineWidth',2);
xlabel('t')
ylim([0,400])
ylabel('F(t)')
title('CD');
xlim([0,tf])
grid on

colors = jet(5);
subplot(2,4,5)
plot(tforwardHFD_treat,out_big_HFD(:,1,1),'color',colors(5,:),'LineWidth',2)
hold on
plot(tforwardHFD_treat,out_big_HFD(:,1,2),'c-','LineWidth',4)
hold on
plot(tforwardHFD_treat,out_big_HFD(:,1,3),'-.','color',colors(3,:),'LineWidth',5);
hold on
plot(tforwardHFD_treat,out_big_HFD(:,1,4),'m--','LineWidth',4);
hold on
plot(tforwardHFD_treat,out_big_HFD(:,1,5),':','color',colors(1,:),'LineWidth',4);
hold on
plot(Tx_uncon_treat_HFD,X_uncon_treat_HFD(:,1),'color',colors(5,:),'LineWidth',2);
xlabel('t')
ylim([0,2100])
ylabel('S(t)')
title('HFD');
xlim([0,tf])
grid on

colors = jet(5);
subplot(2,4,6)
plot(tforwardHFD_treat,out_big_HFD(:,2,1),'color',colors(5,:),'LineWidth',2)
hold on
plot(tforwardHFD_treat,out_big_HFD(:,2,2),'c-','LineWidth',4)
hold on
plot(tforwardHFD_treat,out_big_HFD(:,2,3),'-.','color',colors(3,:),'LineWidth',5);
hold on
plot(tforwardHFD_treat,out_big_HFD(:,2,4),'m--','LineWidth',4);
hold on
plot(tforwardHFD_treat,out_big_HFD(:,2,5),':','color',colors(1,:),'LineWidth',4);
hold on
plot(Tx_uncon_treat_HFD,X_uncon_treat_HFD(:,2),'color',colors(5,:),'LineWidth',2);
xlabel('t')
ylim([0,2100])
ylabel('R(t)')
title('HFD');
xlim([0,tf])
grid on

colors = jet(5);
subplot(2,4,7)
plot(tforwardHFD_treat,out_big_HFD(:,3,1),'color',colors(5,:),'LineWidth',2)
hold on
plot(tforwardHFD_treat,out_big_HFD(:,3,2),'c-','LineWidth',4)
hold on
plot(tforwardHFD_treat,out_big_HFD(:,3,3),'-.','color',colors(3,:),'LineWidth',5);
hold on
plot(tforwardHFD_treat,out_big_HFD(:,3,4),'m--','LineWidth',4);
hold on
plot(tforwardHFD_treat,out_big_HFD(:,3,5),':','color',colors(1,:),'LineWidth',4);
hold on
plot(Tx_uncon_treat_HFD,X_uncon_treat_HFD(:,3),'color',colors(5,:),'LineWidth',2);
xlabel('t')
ylim([0,1500])
ylabel('E(t)')
title('HFD');
xlim([0,tf])
grid on
legend(strcat('p_{HFD}= ',num2str(pHFD_test')),'Location','NorthEast','Orientation','vertical')

subplot(2,4,8)
plot(tforwardHFD_treat,out_big_HFD(:,4,1),'color',colors(5,:),'LineWidth',2)
hold on
plot(tforwardHFD_treat,out_big_HFD(:,4,2),'c-','LineWidth',4)
hold on
plot(tforwardHFD_treat,out_big_HFD(:,4,3),'-.','color',colors(3,:),'LineWidth',5);
hold on
plot(tforwardHFD_treat,out_big_HFD(:,4,4),'m--','LineWidth',4);
hold on
plot(tforwardHFD_treat,out_big_HFD(:,4,5),':','color',colors(1,:),'LineWidth',4);
hold on
plot(Tx_uncon_treat_HFD,X_uncon_treat_HFD(:,4),'color',colors(5,:),'LineWidth',2);
xlabel('t')
ylim([0,400])
ylabel('F(t)')
title('HFD');
xlim([0,tf])
grid on

end

function dy = Model_CD(y,params) % DE

global mu K k3CD m m2CD a1 n alpha k1 k2CD pCD r a2 l a3 cc

dy = zeros(4,1);

SCD = y(1);
RCD = y(2);
ECD = y(3);
FCD = y(4);

dy(1) = k1*SCD*(1- ((SCD+m*RCD)/K))*((ECD^n)/(a1+(ECD^n)))- cc*((a3^l)/(a3^l + ECD^l))*SCD - cc*((a2^l)/(a2^l + ECD^l))*SCD;
dy(2) = k3CD*RCD*(1- ((SCD+m*RCD)/K)) + cc*((a3^l)/(a3^l + ECD^l))*SCD ;
dy(3) = pCD*r*FCD - mu*ECD ;
dy(4) = k2CD*FCD*(1- (FCD*m2CD))- alpha*(SCD+RCD)*FCD;

end

function dy = Model_HFD(y,params) % DE

global mu K k3HFD m m2HFD a1 n alpha k1 k2HFD pHFD r a2 l a3 cc

dy = zeros(4,1);

SHFD = y(1);
RHFD = y(2);
EHFD = y(3);
FHFD = y(4);

dy(1) = k1*SHFD*(1- ((SHFD+m*RHFD)/K))*((EHFD^n)/(a1+(EHFD^n)))- cc*((a3^l)/(a3^l + EHFD^l))*SHFD- cc*((a2^l)/(a2^l + EHFD^l))*SHFD;
dy(2) = k3HFD*RHFD*(1- ((SHFD+m*RHFD)/K)) + cc*((a3^l)/(a3^l + EHFD^l))*SHFD ;
dy(3) = pHFD*r*FHFD - mu*EHFD;
dy(4) = k2HFD*FHFD*(1- (FHFD*m2HFD))- alpha*(SHFD+RHFD)*FHFD;
end
%%
