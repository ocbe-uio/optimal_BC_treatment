%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to simulate Fig.2
%
% Author: Tugba Akman Date: Jan 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data are taken from
% Obesity-Activated Adipose-Derived Stromal Cells Promote
% Breast Cancer Growth and Invasion
% Neoplasia (2018) 20, 1161–1174
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main_Fig2

clear
close all
clc
format long

set(0, 'defaultaxesfontsize',12)


global tumor_CD_tdata tumor_CD_ydata tumor_HFD_tdata tumor_HFD_ydata NumbOfASC_CD_ydata NumbOfASC_HFD_ydata...
    diam_CD diam_HFD Fat_CD_ydata Fat_CD_tdata Fat_HFD_ydata Fat_HFD_tdata tforward tmeasure_Fat tmeasure_tumor...
    initial_cond E0_CD E0_HFD T0 F0_CD F0_HFD weight_tumor_CD weight_tumor_HFD weight_fat_CD weight_fat_HFD...
    mu m1 n k1 a1 r alpha


%Data
Tumor_data_matrix_CD = [65.4	113	179.5	267.9	179.5	113	179.5	113	65.4	113;
    381.5	267.9	381.5	267.9	381.5	267.9	381.5	381.5	267.9	267.9;
    523.3	381.5	523.3	381.5	523.3	381.5	523.3	523.3	523.3	523.3];

Tumor_data_matrix_HFD = [381.5	381.5	381.5	381.5	267.9	179.5	381.5	267.9	381.5	381.5	267.9;
    904.3	904.3	1149.8	696.6	696.6	381.5	696.6	381.5	696.6	523.3	696.6;
    1149.8	1149.8	1436	904.3	904.3	523.3	1436	904.3	1766.3	1149.8	1149.8];

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

%statistical tools
[Fat_data_matrix_CD,Fat_data_matrix_HFD]=Computation_FatData;
[stdev1_CD,weight_tumor_CD,stdev1_fat_CD,weight_fat_CD] = statistical_tools_CD(Tumor_data_matrix_CD,Fat_data_matrix_CD);
[stdev2_HFD,weight_tumor_HFD,stdev2_fat_HFD,weight_fat_HFD] = statistical_tools_HFD(Tumor_data_matrix_HFD,Fat_data_matrix_HFD);

% Time discretization
tmeasure_Fat = 1501; % the points in the solution of F corresponding to the t values of tdata
tmeasure_tumor = [1001,1301,1501]'; % the points in the solution of T corresponding to the t values of tdata

% Parameters in the model
mu = 5.94;%
m1 = 1/2000;
n = 1;
k1 = 0.586967;
a1 = 59.0927;
E0_CD = 175.143;
E0_HFD = 1293.98;
r = 20.8391;
alpha = 2.21427e-05;
F0_CD = 49.923;
F0_HFD = 368.820;
params = [];

%%
% solve the system on a larger interval and plot the solution
tforward = 0:0.01:20; % t mesh for the solution of the DE
initial_cond = [T0 E0_CD F0_CD T0 E0_HFD F0_HFD];
[~, y_r] = ode15s(@(t,y)Model_combined(y,params),tforward,initial_cond);

%%
colormap lines

figure(33)
yyaxis left
f1=plot(tforward,y_r(:,1),'LineWidth',2);
hold on;
errorbar(tumor_CD_tdata,tumor_CD_ydata,stdev1_CD,'ro','LineWidth',2);
hold on
f3=plot(tumor_CD_tdata,tumor_CD_ydata,'ro','LineWidth',2);
hold on
f4=plot(tforward,y_r(:,3),'LineWidth',2);
hold on;
errorbar(Fat_CD_tdata,Fat_CD_ydata,stdev1_fat_CD,'kx--','LineWidth',2);
hold on
f5=plot(Fat_CD_tdata,Fat_CD_ydata,'kx','LineWidth',2);
ylabel('T(t), F(t)')
axis([0 20 0 2000])
yyaxis right
f6=plot(tforward,y_r(:,2),':','LineWidth',2);
ylabel('E(t)')
xlabel('t (days)')
axis([0 20 150 1300])
title('Simulation results for CD')
grid on
legend([f1,  f3,  f4,f5,f6],{'T','Data_T','F','Data_F','E'},'Location','Best')


figure(34)
yyaxis left
f1=plot(tforward,y_r(:,4),'LineWidth',2);
hold on;
errorbar(tumor_HFD_tdata,tumor_HFD_ydata,stdev2_HFD,'ro','LineWidth',2);
hold on
f3=plot(tumor_HFD_tdata,tumor_HFD_ydata,'ro','LineWidth',2);
hold on
f4=plot(tforward,y_r(:,6),'LineWidth',2);
hold on;
errorbar(Fat_HFD_tdata,Fat_HFD_ydata,stdev2_fat_HFD,'kx--','LineWidth',2);
hold on
f5=plot(Fat_HFD_tdata,Fat_HFD_ydata,'kx','LineWidth',2);
ylabel('T(t), F(t)')
axis([0 20 0 2000])
yyaxis right
f6=plot(tforward,y_r(:,5),':','LineWidth',2);
ylabel('E(t)')
xlabel('t (days)')
axis([0 20 150 1300])
title('Simulation results for HFD')
grid on
legend([f1,  f3,  f4,f5,f6],{'T','Data_T','F','Data_F','E'},'Location','Best')

end

%%
function dy = Model_combined(y,params) % DE

global mu m1 k1 a1 n r alpha

dy = zeros(6,1);


TCD = y(1);
ECD = y(2);
FCD = y(3);
THFD = y(4);
EHFD = y(5);
FHFD = y(6);

dy(1) = k1*TCD*(1- (TCD*m1))*((ECD^n)/(a1+(ECD^n)));
dy(2) = r*FCD - mu*ECD ;
dy(3) = - alpha*TCD*FCD;
dy(4) = k1*THFD*(1- (THFD*m1))*((EHFD^n)/(a1+(EHFD^n)));
dy(5) = r*FHFD - mu*EHFD;
dy(6) = - alpha*THFD*FHFD;

end



