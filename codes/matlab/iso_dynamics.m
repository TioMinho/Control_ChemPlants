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
%cd /home/minhotmog/Dropbox/Research/TCC/Codes/
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
        %%
        217,83,79     % 6 - Bootstrap Red
        91,192,222    % 7 - Bootstrap Light Blue
        92,184,92     % 8 - Bootstrap Green
        66,139,202    % 9 - Bootstrap Blue
        255,167,0     % 10 - Google Yellow    
       ]/255;  

%% %%%%%%%%%%%%
%  ISOTHERMAL CSTR 
%  %%%%%%%%%%%%%%
%% Model Loading %%
 run isothermal_cstr/iso_model.m
% load('data/iso_cstr_model.mat')

%% Model Non-Linear Simulation %%
% - Simulation Parameters
% Time
t = (0:0.01:3.99)';                                         
% Initial Conditions
X_0 = [0 0 0 0];
% Input Signal
U = ones(1, numel(t))*3.03;
 
% - Simulation of the Outputs
[~, y] = simulate(iso_cstr.sysVar, t, U, X_0);
    
% - Visualization of the Simulation
figure(1);
subplot(1,2,1), plot(t, U, 'linewidth', 1.5), title("Input Signal")
xlabel("Time (min)"), ylabel("Input Flow-rate (m^3/min)")
ylim([0, 5])
grid()

subplot(1,2,2), p = plot(t, y); title("Output Signals")
p(1).LineWidth = 1.5; p(2).LineWidth = 1.5;
p(3).LineStyle='--'; p(4).LineStyle='--';
legend("C_A", "C_B", "C_C", "C_D") 
xlabel("Time (min)"), ylabel("Outflow Concentration (mol/l)")
grid()

% - Exporting the Visualization to an Image
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('-bestfit', 'isothermal_cstr/simulation/isoCSTR_sim_01', '-dpdf', '-r300')

%% Model Linearized Simulation %%
% - Simulation Parameters
% Linear Model Index
idx = 25;
% Time
t = (0:0.01:1.99)';                                         
% Initial Conditions
U_0 = iso_cstr.oper.U(idx,:); X_0 = iso_cstr.oper.X(idx,:);                                     
% Input Signal
U = [ones(1, numel(t)/4)*1  ones(1, numel(t)/4)*1  ones(1, numel(t)/4)*1 ones(1, numel(t)/4)*1];
 
% Linear Model
A = iso_cstr.ss_model.A_full(X_e, 3.03);     B = iso_cstr.ss_model.B_full(X_e, 3.03);
C = eye(4);          D = zeros(4,1);
iso_cstr_lin = ss(A, B, C, D);

% - Simulation of the Outputs
[~, y] = simulate(iso_cstr.sysVar, t, 3.03+U, [0 0 0 0]);
[y_lin, ~, ~] = lsim(iso_cstr_lin, U, t, [10 10 10 10]/3);
[y_lin_nat, ~, ~] = lsim(iso_cstr_lin, zeros(1,numel(t)), t, [10 10 10 10]/3);
y_lin_forc = y_lin - y_lin_nat;

% - Visualization of the Simulation
figure(2);
subplot(1,2,1), plot(t, U, 'k-', 'linewidth', 1.5), title("Input Signal")
xlabel("Time (min)"), ylabel("u (m^3/min)")
ylim([0, 5])

% subplot(4,2,2), p = plot(t, y(1,:), 'linewidth', 1.5, 'color', cpal(6,:)); hold on; ylabel("x_1 (mol/l)"); title("State Response")
% subplot(4,2,4), p = plot(t, y(2,:), 'linewidth', 1.5, 'color', cpal(8,:)); hold on; ylabel("x_2 (mol/l)")
% subplot(4,2,6), p = plot(t, y(3,:), 'linewidth', 1.5, 'color', cpal(9,:)); hold on; ylabel("x_3 (mol/l)")
% subplot(4,2,8), p = plot(t, y(4,:), 'linewidth', 1.5, 'color', cpal(10,:)); hold on; ylabel("x_4 (mol/l)")
xlabel("Time (min)"), 

subplot(1,2,2), plot(t, y_lin(:,1),'linewidth', 1.5, 'linestyle', '-', 'color', cpal(6,:)); hold on; ylabel("x_1 (mol/l)"); title("State Response")
subplot(1,2,2), plot(t, y_lin_nat(:,1),'linewidth', 1.5, 'linestyle', '--', 'color', cpal(6,:)); hold on
subplot(1,2,2), plot(t, y_lin_forc(:,1),'linewidth', 1.5, 'linestyle', ':', 'color', cpal(6,:)); 
% subplot(4,2,4), plot(t, y_lin(:,2)+X_e(2),'linewidth', 1.5, 'linestyle', '--', 'color', cpal(8,:)); hold on
% subplot(4,2,6), plot(t, y_lin(:,3)+X_e(3),'linewidth', 1.5, 'linestyle', '--', 'color', cpal(9,:)); hold on
% subplot(4,2,8), plot(t, y_lin(:,4)+X_e(4),'linewidth', 1.5, 'linestyle', '--', 'color', cpal(10,:))
xlabel("Time (min)"), 

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'isothermal_cstr/simulation/isoCSTR_sim_lin_99', '-dpdf', '-r300')

%% Operation Points Visualization %%
figure(3); 
hold on
for n = 1:iso_cstr.oper.size-1
    plot(iso_cstr.oper.U([n n+1]), iso_cstr.oper.X([n n+1],2), 'color', ccmap.linear(n,:), 'linewidth', 2);
end
title("Stationary Points - C_B")
xlabel("Input Flow-rate (m^3/min)"), ylabel("Outflow Concentration (mol/m^3)")
grid()

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'isothermal_cstr/simulation/isoCSTR_sim_op_01', '-dpdf', '-r300')

%% System Modes Visualization %%
t = 0:0.1:3;
figure(4);
for i = 1:iso_cstr.oper.size
    tic
    subplot(1,2,1), plot(t, subs(iso_cstr.modes(i,1), t), 'color', ccmap.linear(i,:))
    title("First Mode"), xlabel("Time (min)"), ylabel("e^{\lambda_1 t}"), hold on
    
    subplot(1,2,2), plot(t, subs(iso_cstr.modes(i,2), t), 'color', ccmap.linear(i,:))
    title("Second Mode"), xlabel("Time (min)"), ylabel("e^{\lambda_2 t}"), hold on
    
    drawnow
    toc
end

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'isothermal_cstr/simulation/isoCSTR_sim_modes_01', '-dpdf', '-r300')

%% Transfer Matrices Visualization %%
t = 0:0.1:4;
figure(5);
for i = 1:iso_cstr.oper.size
    tic
    subplot(2,2,1), plot(t, subs(iso_cstr.trf_matrix(1,1,i), t), 'color', ccmap.linear(i,:)), hold on
    subplot(2,2,2), plot(t, subs(iso_cstr.trf_matrix(1,2,i), t), 'color', ccmap.linear(i,:)), hold on
    subplot(2,2,3), plot(t, subs(iso_cstr.trf_matrix(2,1,i), t), 'color', ccmap.linear(i,:)), hold on
    subplot(2,2,4), plot(t, subs(iso_cstr.trf_matrix(2,2,i), t), 'color', ccmap.linear(i,:)), hold on
    drawnow
    toc
end

subplot(2,2,1), title("X_1"), ylabel("X_1")
subplot(2,2,2), title("X_2")
subplot(2,2,3), ylabel("X_2"), xlabel("Time (min)")
subplot(2,2,4), xlabel("Time (min)")

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'isothermal_cstr/simulation/isoCSTR_sim_expAt_01', '-dpdf', '-r300')

%% Stability Analysis - Pole-Zero Mapping %%
figure(6);
for i = 1:iso_cstr.oper.size
    tic
    plot(iso_cstr.poles(i,:), 0, 'x', 'color', ccmap.linear(i,:), 'markersize', 10); hold on
    drawnow
    toc
end
title("Pole-Zero Map")
xlabel("Real Axis (seconds^{-1})")
ylabel("Imaginary Axis (seconds^{-1})")
sgrid()

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'isothermal_cstr/simulation/isoCSTR_stab_pz_01', '-dpdf', '-r300')

%% Stability Analysis - Bode Plots %%
figure(7);
for i = 1:iso_cstr.oper.size
    tic
    bodeplot(ss(iso_cstr.ss_model.A(i), iso_cstr.ss_model.B(i), iso_cstr.ss_model.C, iso_cstr.ss_model.D), 'b'), hold on;

    lineHandle = findobj(gcf,'Type','line','-and','Color','b');
    set(lineHandle,'Color',ccmap.linear(i,:));
    
    drawnow
    toc
end
grid()

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'isothermal_cstr/simulation/isoCSTR_stab_bode_01', '-dpdf', '-r300')

%% Stability Analysis - Nyquist Plots %%
figure(8);
for i = 1:iso_cstr.oper.size
    tic
    nyquistplot(ss(iso_cstr.ss_model.A(i), iso_cstr.ss_model.B(i), iso_cstr.ss_model.C, iso_cstr.ss_model.D), 'b'), hold on;

    lineHandle = findobj(gcf,'Type','line','-and','Color','b');
    set(lineHandle,'Color',ccmap.linear(i,:));
    
    drawnow
    toc
end
grid()

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'isothermal_cstr/simulation/isoCSTR_stab_nyquist_01', '-dpdf', '-r300')
