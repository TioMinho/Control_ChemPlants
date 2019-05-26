%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamical Modelling and Intelligent Control of Chemical Reactive Systems
% 
% Author: Otacilio Bezerra Leite Neto
%
%   exo_dynamics.m
%       This script contains the instructions for running and visualize experiments 
%       over the dynamics of the Exothermal Continuous-Stirred Tank (CSTR) system.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preamble %%
clc; close all; clear all;
% cd /home/minho/Documents/Minho/Control_ChemPlants/
% cd /home/minhotmog/Dropbox/Research/TCC/Codes/

% Sets the Default Rendere to tbe the Painters
set(0, 'DefaultFigureRenderer', 'painters');

% Some colors
load('data/ccmap.mat');
cpal = [209 17 65;    % 1 - Metro Red
        0 177 89;     % 2 - Metro Green
        0 174 219;    % 3 - Metro Blue
        243 119 53;   % 4 - Metro Orange
        255 196 37;   % 5 - Metro Yellow
        %%
        217,83,79;     % 6 - Bootstrap Red
        91,192,222;    % 7 - Bootstrap Light Blue
        92,184,92;     % 8 - Bootstrap Green
        66,139,202;    % 9 - Bootstrap Blue
        255,167,0;     % 10 - Google Yellow
        lines(10)*255
       ]/255;  
                          
%% %%%%%%%%%%%%
%  EXOTHERMAL CSTR %%
%  %%%%%%%%%%%%%%
%% Model Loading %%
run exothermal_cstr/exo_model.m
%load('data/exo_cstr_model.mat')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%
%  Model Non-Linear Simulation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Parameters
% Time
t = (0:0.01:1)'; T = numel(t);
% Initial Conditions
X_0 = [1.235 0.9 134.14 128.95 2.4192 0.2730];
% Input Signal
U = [ones(1,floor(T/3))*18.83 ones(1,T-2*floor(T/3))*30 ones(1,floor(T/3))*18.83
    ones(1,floor(T/3))*-4495 ones(1,T-2*floor(T/3))*-6495 ones(1,floor(T/3))*-4495];   
 
% Simulation of the Outputs
[~, y] = simulate(exo_cstr.sysVar, t, U, X_0);
    
%% Visualization
figure(1);
subplot(2,2,1), plot(t, U(1,:), 'linewidth', 1.5, 'color', [0.3 0.3 0.3])
xlabel("Time (min)"), ylabel("Flow-rate (1/hr)")

subplot(2,2,2)
p = plot(t, y(1, :), 'linewidth', 1.5, 'color', cpal(11,:)); hold on
p = plot(t, y(2, :), 'linewidth', 1.5, 'color', cpal(12,:)); hold on
p = plot(t, y(5, :), 'linestyle', '--', 'linewidth', 1, 'color', cpal(13,:)); hold on
p = plot(t, y(6, :), 'linestyle', '--', 'linewidth', 1, 'color', cpal(14,:));
xlabel("Time (min)"), ylabel("Concentration (mol/l)")

subplot(2,2,3), plot(t, U(2,:), 'linewidth', 1.5, 'color', [0.3 0.3 0.3])
xlabel("Time (min)"), ylabel("Cooling Capacity (kJ/hr)")

subplot(2,2,4), 
plot(t, y(3, :), 'linewidth', 1.5, 'color', cpal(17,:)); hold on
plot(t, y(4, :), 'linewidth', 1.5, 'color', cpal(16,:));
xlabel("Time (min)"), ylabel("Temperatures (ºC)")

% - Exporting the Visualization to an Image
timeNow = datetime('now', 'TimeZone', 'local', 'Format', 'dMMMy_HHmmssZ');
figname = "report_experiments/figs/exoSim_dynamics_" + char(timeNow);
fig = gcf; fig.PaperPositionMode = 'auto'; 
print('-bestfit', figname, '-dpdf', '-r300')
system("pdfcrop " + figname + ".pdf " + figname + ".pdf");

%% %%%%%%%%%%%%%%%%%%%%%%%%%%
%  Model Linearized Simulation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Simulation Parameters
% Time
t = (0:0.1:10)'; T = numel(t);

% Initial Conditions
U_ss = [18.83 -4495]; X_ss = [1.235 0.9 134.14 128.95];

% Input Signal
U = [ones(1,floor(T/3))*U_ss(1) ones(1,T-2*floor(T/3))*22 ones(1,floor(T/3))*U_ss(1);
     ones(1,floor(T/3))*U_ss(2) ones(1,T-2*floor(T/3))*-5000 ones(1,floor(T/3))*U_ss(2)];

U = U + randn(2,T) .* (0.01*U);
 
% Linear Model
A = exo_cstr.ss_model.A_l(X_ss, U_ss); B = exo_cstr.ss_model.B_l(X_ss, U_ss);
C = exo_cstr.ss_model.C;               D = exo_cstr.ss_model.D;

exo_cstr_lin = ss(A, B, C, D);

% - Simulation of the Outputs
[~, y] = simulate(exo_cstr.sysVar, t, U, [X_ss 2.4192 0.2730]);
[y_lin, ~, ~] = lsim(exo_cstr_lin, U-U_ss', t);

y_lin = y_lin+X_ss;

%% Visualization
figure(2);
subplot(2,2,1), plot(t, U(1,:), 'linewidth', 1.5, 'color', [0.3 0.3 0.3])
xlabel("Time (min)"), ylabel("Flow-rate (1/hr)")

subplot(2,2,2)
p = plot(t, y(1, :), 'linewidth', 1, 'color', cpal(11,:)); hold on
p = plot(t, y(2, :), 'linewidth', 1, 'color', cpal(12,:)); hold on
p = plot(t, y_lin(:, 1), 'linestyle', '--', 'linewidth', 1.5, 'color', cpal(11,:)); hold on
p = plot(t, y_lin(:, 2), 'linestyle', '--', 'linewidth', 1.5, 'color', cpal(12,:));
xlabel("Time (min)"), ylabel("Concentration (mol/l)")

subplot(2,2,3), plot(t, U(2,:), 'linewidth', 1.5, 'color', [0.3 0.3 0.3])
xlabel("Time (min)"), ylabel("Cooling Capacity (kJ/hr)")

subplot(2,2,4), 
plot(t, y(3, :), 'linewidth', 1, 'color', cpal(17,:)); hold on
plot(t, y(4, :), 'linewidth', 1, 'color', cpal(16,:)); hold on
plot(t, y_lin(:, 3), 'linestyle', '--', 'linewidth', 1.5, 'color', cpal(17,:)); hold on
plot(t, y_lin(:, 4), 'linestyle', '--', 'linewidth', 1.5, 'color', cpal(16,:));
xlabel("Time (min)"), ylabel("Temperatures (ºC)")

% - Exporting the Visualization to an Image
timeNow = datetime('now', 'TimeZone', 'local', 'Format', 'dMMMy_HHmmssZ');
figname = "report_experiments/figs/exoSim_dynamics_" + char(timeNow);
fig = gcf; fig.PaperPositionMode = 'auto'; 
print('-bestfit', figname, '-dpdf', '-r300')
system("pdfcrop " + figname + ".pdf " + figname + ".pdf");

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Time Response Analysis (Modes) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the system modes in a symbolic expression
syms t
[V, J] = jordan(A);
e_At = V*expm(J*t)*pinv(V);

% Simulate the contribution of each mode
t = 0:0.001:1.5;

%% Visualization
figure(4);

for j = 1:4 
for k = 1:4
    tic
    subplot(4,4,j+4*(k-1)), plot(t, subs(e_At(j,k), t), 'color', cpal(11,:))
    toc
end
end

% - Exporting the Visualization to an Image
timeNow = datetime('now', 'TimeZone', 'local', 'Format', 'dMMMy_HHmmssZ');
figname = "report_experiments/figs/exoSim_dynamics_" + char(timeNow);
fig = gcf; fig.PaperPositionMode = 'auto'; 
print('-bestfit', figname, '-dpdf', '-r300')
system("pdfcrop " + figname + ".pdf " + figname + ".pdf");