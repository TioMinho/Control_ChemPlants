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
% - Simulation Parameters
% Linear Model Index
idx = 25;
% Time
t = (0:0.05:2)';                                         
% Initial Conditions
U_0 = iso_cstr.oper.U(idx,:); X_0 = iso_cstr.oper.X(idx,:);
 
% Linear Model
A = iso_cstr.ss_model.A(idx); B = iso_cstr.ss_model.B(idx);
C = iso_cstr.ss_model.C;        D = iso_cstr.ss_model.D;
K = lqr(A, B, [2 0; 0 2],  [2]);
iso_cstr_lin = ss(A-K*B, B, iso_cstr.ss_model.C, iso_cstr.ss_model.D);

% - Simulation of the Outputs
[y_lin, ~, U, ~] = lqr_sim(A, B, C, D, [3 0; 0 7] , [5], [15, 0] - X_0, t);  
[~, y] = simulate(iso_cstr.model, t, U_0+U, [15 0 0 0]);

% [y_lin2, ~] = lsim((A-B*K), B, C, D, zeros(1,numel(t)), t, [10, 2] - X_0);
    
% - Visualization of the Simulation
figure(1);
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

subplot(2,2,4), plot(t, y_lin+X_0, 'linewidth', 1.5);
line([0, t(end)], [X_0(1), X_0(1)], 'linestyle', '--', 'color', 'black'); 
line([0, t(end)], [X_0(2), X_0(2)], 'linestyle', '--', 'color', 'black')
title("Linearized CSTR"), xlabel("Time (min)"), ylabel("Outflow Concentration (mol/l)")
legend("C_A", "C_B") 
grid()

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'isothermal_cstr/simulation/isoCSTR_lqr-ih_04', '-dpdf', '-r300')

%% Discrete-Time Finite Horizon Linear Quadratic Regulator %%
% - Simulation Parameters
% Linear Model Index
idx = 25;
figure(2);
for idx = 25:25
    tic
    % Time
    t = (0:0.05:2)';                                         
    % Initial Conditions
    U_0 = iso_cstr.oper.U(idx,:); X_0 = iso_cstr.oper.X(idx,:);

    % Linear Model
    A = iso_cstr.ss_model.A(idx); B = iso_cstr.ss_model.B(idx);
    C = iso_cstr.ss_model.C;        D = iso_cstr.ss_model.D;
    iso_cstr_lin = ss(A-K*B, B, iso_cstr.ss_model.C, iso_cstr.ss_model.D);

    % - Simulation of the Outputs
    [y_lin, ~, U, ~] = lqr_fh_sim(A, B, C, D, [10 0; 0 40] , [5], 2, [10, 2] - X_0, t);  
    [~, y] = simulate(iso_cstr.model, t, U_0+U, [10, 2 0 0]);

    % - Visualization of the Simulation

    subplot(2,2,1), plot(t, U_0+U, 'linewidth', 1.5); 
    title("Input Signal")
    xlabel("Time (min)"), ylabel("Input Flow-rate (m^3/min)")
    grid()

    subplot(2,2,2), p = plot(t, y); 
    title("Non-Linear CSTR")
    p(1).LineWidth = 1.5; p(2).LineWidth = 1.5;
    p(3).LineStyle='--'; p(4).LineStyle='--';
    legend("C_A", "C_B", "C_C", "C_D") 
    xlabel("Time (min)"), ylabel("Outflow Concentration (mol/l)")
    grid()

    subplot(2,2,4), plot(t, y_lin,'linewidth', 1.5); 
    title("Linearized CSTR")
    legend("C_A", "C_B") 
    xlabel("Time (min)"), ylabel("Outflow Concentration (mol/l)")
    grid()

    drawnow
    toc
end

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'isothermal_cstr/simulation/isoCSTR_lqr-dh_01', '-dpdf', '-r300')