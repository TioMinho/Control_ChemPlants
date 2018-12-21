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

%% Continuous-Time Infinite Horizon Linear Quadratic Regulator %%
clc
% - Simulation Parameters
% Linear Model Index
idx = 25;
% Time
t = (0:0.05:19.95)';
% Initial Conditions
U_0 = iso_cstr.oper.U(idx,:); X_0 = iso_cstr.oper.X(idx,:);

% Reference Signall
r = [ones(1,50)*X_0(1) ones(1,150)*iso_cstr.oper.X(2,1) ones(1,150)*iso_cstr.oper.X(50,1) ones(1,50)*iso_cstr.oper.X(idx,1);
       ones(1,50)*X_0(2) ones(1,150)*iso_cstr.oper.X(2,2) ones(1,150)*iso_cstr.oper.X(50,2) ones(1,50)*iso_cstr.oper.X(idx,2);];

% Disturbance signal
w = randn(numel(t), 1) * 0.1;

% Linear Model
A = iso_cstr.ss_model.A(idx);   B = iso_cstr.ss_model.B(idx);
C = iso_cstr.ss_model.C;        D = iso_cstr.ss_model.D;

% Controller and Observer
Q = diag([20, 30]);
R = [10];
L = [0 0; 
	   0 0];

% - Simulation of the Outputs
[~, yout, xout, uout] = simulate(iso_cstr, idx, t, r, X_0, 'lqg', Q, R, numel(t), L, w);
    
% - Visualization of the Simulation
figure(1);
subplot(1,2,1), plot(t, U_0+uout, 'linewidth', 1.5), title("Input Signal")
xlabel("Time (min)"), ylabel("Input Flow-rate (m^3/min)")
grid()

subplot(1,2,2)
plot(t, r, 'linestyle', '--', 'color', 'black'); hold on;

p = plot(t, yout+X_0', 'linewidth', 1.5, 'linestyle', '--');
set(p, {'color'}, {cpal(3,:); cpal(4,:)});

p = plot(t, xout, 'linewidth', 1.5); hold on
set(p, {'color'}, {cpal(3,:); cpal(4,:)});

title("Isothermal CSTR"), xlabel("Time (min)"), ylabel("Outflow Concentration (mol/l)")
legend("C_{Ar}", "C_{Ar}", "C_{Ae}", "C_{Be}", "C_{A}", "C_{B}") 
grid()

% - Exporting the Visualization to an Image
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('-bestfit', 'isothermal_cstr/simulation/isoCSTR_lqr-ih_04', '-dpdf', '-r300')

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