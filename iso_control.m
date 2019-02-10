%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamical Modelling and Intelligent Control of Chemical Reactive Systems
% 
% Author: Otacilio Bezerra Leite Neto
%
%   iso_dynamics.m
%       This script contains the instructions for running and visualize experiments 
%       over the dynamics of the Isothermal Continuous-Stirred Tank (CSTR) system.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preamble %%
% cd /home/omenezes/Documents/Minho/Control_ChemPlants/
cd /home/minhotmog/Dropbox/Research/TCC/Codes/
clc; clear all; close all;

% Sets the Default Rendere to tbe the Painters
set(0, 'DefaultFigureRenderer', 'painters');

% Some colors
load('data/ccmap.mat');
cpal = [209 17 65;    % 1 - Metro Red
             0 177 89;     % 2 - Metro Green
            0 174 219;    % 3 - Metro Blue
           243 119 53;   % 4 - Metro Orange
           255 196 37;   % 5 - Metro Yellow
                              ]/255;  

%% %%%%%%%%%%%%
%  ISOTHERMAL CSTR 
%  %%%%%%%%%%%%%%
%% Model Loading %%
run isothermal_cstr/iso_model.m
% load('data/iso_cstr_model.mat')

%% Continuous-Time Infinite Horizon Linear Quadratic Regulator %%
% - Simulation Parameters
% Linear Model Index
idx = 25;
% Time
t = (0:0.05:20.95)';
% Initial Conditions
U_0 = iso_cstr.oper.U(idx,:); X_0 = iso_cstr.oper.X(idx,:);

% Reference Signal
r = [%ones(1,50)*iso_cstr.oper.X(idx,1) ones(1,160)*iso_cstr.oper.X(15,1) ones(1,160)*iso_cstr.oper.X(35,1) ones(1,50)*iso_cstr.oper.X(idx,1);
       ones(1,50)*iso_cstr.oper.X(idx,2) ones(1,160)*iso_cstr.oper.X(15,2) ones(1,160)*iso_cstr.oper.X(35,2) ones(1,50)*iso_cstr.oper.X(idx,2)];

% Disturbance signal
w = randn(numel(t), 1) * .1;             % Process Noise
z = randn(numel(t), 2) .* [.075 .05];      % Measurement Noise

% Linear Model
A = iso_cstr.ss_model.A(idx);   B = iso_cstr.ss_model.B(idx);
C = iso_cstr.ss_model.C;          D = iso_cstr.ss_model.D;

iso_cstr.ss_model.C = [0 1];
iso_cstr.ss_model.D = [0];
iso_cstr.sizeY = 1;

% Controller and Observer
Q1 = diag([20, 20, 1e4]);
R1 = diag([5]);

L = [0 0; 0 0];

% - Simulation of the Outputs
[~, yout, xout, uout] = simulate(iso_cstr, idx, t, r, X_0, 'lqgi', Q1, R1, 'inf', L, w, z);
mts1 = metrics(yout, r, t)

% - Visualization of the Simulation
figure(1);
subplot(1,2,1)
plot(t, U_0+uout, 'linewidth', 1.5, 'color', cpal(1,:)); 
xlabel("Time (min)"), ylabel("Input Flow-rate (m^3/min)"),  title("Input Signal")
grid()

subplot(2,2,2)
% plot(t, r(1,:), 'linestyle', '--', 'color', 'black'); hold on;
scatter(t, xout(1,:), 'marker', 'x', 'MarkerEdgeColor', cpal(3,:)); hold on; 
plot(t, xout(1,:)-z(:,1)', 'linewidth', 1.5, 'linestyle', '--', 'color', cpal(1,:));
plot(t, yout(1,:)+X_0(:,1)', 'linewidth', 1.5, 'linestyle', '-', 'color', cpal(2,:)); hold on;
xlabel("Time (min)"), ylabel("Outflow Concentration - C_a (mol/l)"), title("Isothermal CSTR")
legend('y(t)', 'x(t)', 'x\_hat(t)')
grid()

subplot(2,2,4), 
plot(t, r, 'linestyle', '--', 'color', 'black'); hold on;
scatter(t, xout(2,:), 'marker', 'x', 'MarkerEdgeColor', cpal(3,:)); hold on; 
plot(t, xout(2,:)-z(:,2)', 'linewidth', 1.5, 'linestyle', '--', 'color', cpal(1,:));
plot(t, yout(2,:)+X_0(:,2)', 'linewidth', 1.5, 'linestyle', '-', 'color', cpal(2,:)); hold on; 
xlabel("Time (min)"), ylabel("Outflow Concentration - C_b (mol/l)")
legend('Reference', 'y(t)', 'x(t)', 'x\_hat(t)')
grid()

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
set(fig,'PaperOrientation', 'landscape');
print('-bestfit', 'isothermal_cstr/simulation/isoCSTR_lqri_01', '-dpdf', '-r300')

%% Discrete-Time Finite Horizon Linear Quadratic Regulator %%
% - Simulation Parameters
% Linear Model Index
idx = 25;

% Time
t = (0:0.05:2)';                                         
% Initial Conditions
U_0 = iso_cstr.oper.U(idx,:); X_0 = iso_cstr.oper.X(idx,:);

% Linear Model
A = iso_cstr.ss_model.A(idx);   B = iso_cstr.ss_model.B(idx);
C = iso_cstr.ss_model.C;        D = iso_cstr.ss_model.D;

[K, P] = lqr_(A, B, [2 0; 0 2], [2], 10);

% - Simulation of the Outputs
[~, y] = simulate(iso_cstr.model, t, U_0+U, [0.1, 10 0 0]);

% - Visualization of the Simulation
figure(2);
subplot(2,2,1), plot(t, U_0+U, 'linewidth', 1.5), title("Input Signal")
xlabel("Time (min)"), ylabel("Input Flow-rate (m^3/min)")
grid()

subplot(2,2,2), 
p = plot(t, y(:,1:2), 'linewidth', 1.5);
line([0, t(end)], [X_0(1), X_0(1)], 'linestyle', '--', 'color', 'black'); 
line([0, t(end)], [X_0(2), X_0(2)], 'linestyle', '--', 'color', 'black')
title("Non-Linear CSTR"), xlabel("Time (min)"), ylabel("Outflow Concentration (mol/l)")
legend("C_A", "C_B") 
grid()

subplot(2,2,4), 
plot(t, y_lin_hat+X_0, 'linewidth', 1.5); hold on; 
plot(t, y_lin+X_0, 'linestyle', '--');
line([0, t(end)], [X_0(1), X_0(1)], 'linestyle', '--', 'color', 'black'); 
line([0, t(end)], [X_0(2), X_0(2)], 'linestyle', '--', 'color', 'black')
title("Linearized CSTR"), xlabel("Time (min)"), ylabel("Outflow Concentration (mol/l)")
legend("C_A^{est.}", "C_B^{est.}", "C_A^{real}", "C_B^{real}") 
grid()

% - Exporting the Visualization to an Image
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('-bestfit', 'isothermal_cstr/simulation/isoCSTR_lqr-dh_05', '-dpdf', '-r300')