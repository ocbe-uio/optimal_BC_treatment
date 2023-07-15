%% Full model 

close all; clear all; clc;

diary('diary_MiceExp.txt')
diary on

set(0, 'defaultaxesfontsize',16)

% initialize full model
Setup

%% Optimization
arSimu(true,true,true);

%arPlot

%return

n_fit = 2000; % number of starts

if n_fit > 1
    arFitLHS(n_fit); % multistart optimization
else
    arFit;
end

% Print optimization results
arPrint;     % display parameter values
arSimu(true,true,true);

% %Profile likelihood
% arPLEInit
% ple

%%
arPlot

% To save
arSave('ModelExp_Test1')

diary off
