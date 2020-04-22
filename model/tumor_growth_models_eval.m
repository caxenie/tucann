% Tumor Characterization using Unsupervised Learning of Mathematical Relations within Breast Cancer Data by Cristian Axenie and Daria Kurz
% prepare environment
clear;
clc; close all; 
% dataset name
dataset = 'Experiment_dataset_*.mat';
load(dataset);
global T M
% data points and time indices
T = 1:DATASET_LEN_ORIG; T = T';
M = sensory_data_orig.y;
%%  ODE integration is sensitive to initial conditions
% some experiments have a scale 3 order of magnitude larger as others
% adjust the initial condition for these two experiments

% Dataset 1: Rodallec, Anne et al, Tumor growth kinetics of human MDA-MB-231
if experiment_dataset == 1
    minv = 0.2;
else
    minv = 0.02;
end
%% evaluate classical models 
[t1g, y1g] = tumor_growth_model_fit(T, M,'Gompertz', minv);
[t1l, y1l] = tumor_growth_model_fit(T, M,'logistic', minv);
[t1v, y1v] = tumor_growth_model_fit(T, M,'vonBertalanffy', minv);
% evaluate the neural model 
t1neuro = T;
y1neuro = neural_model(T)';
% evalute visually 
figure();
set(gcf, 'color', 'w');hold on; box off;
% dataset points
plot(T,M,'r*'); 
% neural model 
plot(t1neuro, y1neuro);
% classical models
plot(t1g,y1g); 
plot(t1l,y1l);
plot(t1v,y1v);
title('Growth models analysis');
legend('Data points', 'Ours', 'Gompertz','Logistic','vonBertalanffy');
legend('boxoff');
box off;
%% boxplot evaluation
figure; set(gcf,'color', 'w'); box off;
% combine the prediction in unified vector
x = [M;y1neuro;y1g;y1l;y1v];
% create a grouping variable
g1 = repmat({'Data'}, length(M), 1);
g2 = repmat({'Ours'}, length(y1neuro), 1);
g3 = repmat({'Gompertz'}, length(y1g), 1);
g4 = repmat({'Logistic'}, length(y1l), 1);
g5 = repmat({'vonBertalanffy'}, length(y1v), 1);
g=[g1;g2;g3;g4;g5];
boxplot(x, g); box off;
%% Evaluate SSE, RMSE, MAPE

% resample
y1g = interp1(1:length(y1g), y1g, linspace(1,length(y1g),length(M)))';
y1l = interp1(1:length(y1l), y1l, linspace(1,length(y1l),length(M)))';
y1v = interp1(1:length(y1v), y1v, linspace(1,length(y1v),length(M)))';

% params as from [Benzekry et al., 2014c]
alfa = 0.84;
sigma = 0.21;

% locals, model sequence, names and param numbers from [Benzekry et al., 2014c]
models = 1:4;
names = {'Ours'; 'Gompertz'; 'Logistic'; 'Bertalanffy'};
param_num = [0, 2, 2, 3]; % param number for each model 

% SSE
SSEn = model_sse(alfa, sigma, M, y1neuro);
SSEg = model_sse(alfa, sigma, M, y1g);
SSEl = model_sse(alfa, sigma, M, y1l);
SSEv = model_sse(alfa, sigma, M, y1v);

%plot comparatively
figure(); set(gcf, 'color', 'w'); 
plot(models, [SSEn, SSEg, SSEl, SSEv],'k*');box off;
ylabel('SSE');
set(gca,'xtick', models, 'xticklabel', names);

% RMSE
RMSEn = model_rmse(alfa, sigma, param_num(1), M, y1neuro);
RMSEg = model_rmse(alfa, sigma, param_num(2), M, y1g);
RMSEl = model_rmse(alfa, sigma, param_num(3), M, y1l);
RMSEv = model_rmse(alfa, sigma, param_num(4), M, y1v);

%plot comparatively
figure(); set(gcf, 'color', 'w'); 
plot(models, [RMSEn, RMSEg, RMSEl, RMSEv],'k*');
ylabel('RMSE');
set(gca,'xtick', models, 'xticklabel', names);box off;

% sMAPE
sMAPEn = mean(2*abs((M-y1neuro))./(abs(M) + abs(y1neuro)));
sMAPEg = mean(2*abs((M-y1g))./(abs(M) + abs(y1g)));
sMAPEl = mean(2*abs((M-y1g))./(abs(M) + abs(y1g)));
sMAPEv = mean(2*abs((M-y1g))./(abs(M) + abs(y1g)));

%plot comparatively
figure(); set(gcf, 'color', 'w'); 
plot(models, [sMAPEn, sMAPEg, sMAPEl, sMAPEv],'k*');box off;
ylabel('sMAPE');
set(gca,'xtick', models, 'xticklabel', names);