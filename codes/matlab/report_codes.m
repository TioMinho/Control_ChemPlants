%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamical Modelling and Intelligent Control of Chemical Reactive Systems
% 
% Author: Otacilio Bezerra Leite Neto
%
%   report_codes.m
%       This script contains the codes used to generate the images displayed in the report 
%       for the monograph. (Diagrams and Illustrations are not included)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preamble %%
% cd /home/minho/Documents/Minho/Control_ChemPlants/
% cd /home/minhotmog/Dropbox/Research/TCC/Codes/
clc; clear all; close all;

% Sets the Default Renderer to tbe the Painters
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
        92,184,92     % 7 - Bootstrap Green
        66,139,202   % 8 - Bootstrap Blue
        255,167,0     % 9 - Google Yellow    
       ]/255;  

%% %%%%%%%%%%%%
%  ISOTHERMAL CSTR 
%  %%%%%%%%%%%%%%
%% Model Loading %%
run isothermal_cstr/iso_model.m
% load('data/iso_cstr_model.mat')

%% Figure 2.3 %%
% - Simulation of the Nonlinear plat to obtain steady-state values
[~, y] = simulate(iso_cstr.sysVar, (0:.1:4)', ones(1,41)*3.03, [0 0 0 0]);
x_o = y(:,end);

% - Real simulation parameters
t = (0:0.01:5.99)';                                          
X_0 = [0 0 0 0];
U = [ones(1, numel(t)/4)*1  ones(1, numel(t)/4)*3  ones(1, numel(t)/4)*2 ones(1, numel(t)/4)*4];

% - State-space model selection
A = iso_cstr.ss_model.A_full(x_o, 3.03);     B = iso_cstr.ss_model.B_full(x_o, 3.03);
C = eye(4);                                                  D = zeros(4,1);
iso_cstr_lin = ss(A, B, C, D);

% - Simulation of the Outputs
[~, y] = simulate(iso_cstr.sysVar, t, U+3.03, X_0);
[y_lin, ~, ~] = lsim(iso_cstr_lin, U, t, -x_o);

% - Visualization of the Simulation
figure(1);

subplot(2,2,1), plot(t, U, 'k-', 'linewidth', 1.5), title("Input Signal")
xlabel("Time (min)"), ylabel("u (m^3/min)")
ylim([0, 5])


for i = 1:4
    subplot(4,2,2*i), plot(t, y(i,:), 'linewidth', 1.5, 'color', cpal(5+i, :)); hold on;
    subplot(4,2,2*i), plot(t, y_lin(:,i)+x_o(i), 'linewidth', 1.5, 'linestyle', '--', 'color', cpal(5+i, :))
    ylabel(strcat("x_", int2str(i), " (mol/l)"))
end
subplot(4,2,2); title("State Response")
subplot(4,2,8); xlabel("Time (min)")

% - Exporting the Visualization to an Image
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('-bestfit', 'isothermal_cstr/simulation/isoCSTR_sim_01', '-dpdf', '-r300')

%% Figure 2.4 %%
% - Simulation of the Nonlinear plat to obtain steady-state values
[~, y] = simulate(iso_cstr.sysVar, (0:.1:4)', ones(1,41)*3.03, [0 0 0 0]);
x_o = y(1:2,end);

% - Real simulation parameters
t = (0:0.01:1.99)';      
U = ones(1,numel(t));
U_zero = zeros(1,numel(t));

% - State-space model selection
A = iso_cstr.ss_model.A_l(x_o, 3.03);     B = iso_cstr.ss_model.B_l(x_o, 3.03);
C = eye(2);                                               D = zeros(2,1);
iso_cstr_lin = ss(A, B, C, D);

% - Simulation of the Outputs
[y_lin, ~, ~]        = lsim(iso_cstr_lin, U, t, [10 10]/3);
[y_lin_nat, ~, ~] = lsim(iso_cstr_lin, U_zero, t, [10 10]/3);
y_lin_forc = y_lin - y_lin_nat;

% - Visualization of the Simulation
figure(2);

subplot(1,2,1), plot(t, U, 'k-', 'linewidth', 1.5), title("Input Signal")
xlabel("Time (min)"), ylabel("u (m^3/min)")
ylim([0, 5])

subplot(1,2,2)
plot(t, y_lin(:,1), 'linewidth', 1.5, 'color', cpal(8, :)); hold on;
plot(t, y_lin_nat(:,1), 'linewidth', 1.5, 'linestyle', '--', 'color', cpal(8, :)); hold on;
plot(t, y_lin_forc(:,1), 'linewidth', 1.5, 'linestyle', ':', 'color', cpal(8, :)); hold on;

ylabel(strcat("x_1 (mol/l)"))
title("State Response")
xlabel("Time (min)")

% - Exporting the Visualization to an Image
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('-bestfit', 'isothermal_cstr/simulation/isoCSTR_sim_01', '-dpdf', '-r300')

%% Figure 2.5 %%
% - Simulation of the Nonlinear plat to obtain steady-state values
[~, y] = simulate(iso_cstr.sysVar, (0:.1:4)', ones(1,41)*3.03, [0 0 0 0]);
x_o = y(1:2,end);

% - Real simulation parameters
t_s = (0:0.01:1.99)';      
U = ones(1,numel(t));
U_zero = zeros(1,numel(t));

% - State-space model selection
A = iso_cstr.ss_model.A_l(x_o, 3.03);     B = iso_cstr.ss_model.B_l(x_o, 3.03);
C = eye(2);                                               D = zeros(2,1);
iso_cstr_lin = ss(A, B, C, D);

% - Calculation of the Matrix Exponential
syms t
e_At = expm(A * t);

% - Visualization of the Matrix Exponential
figure(3);
for i = 1:2
    for j = 1:2
        subplot(2,2,(i-1)*2+j)
        plot(t_s, double(subs(e_At(i,j), t_s)), 'linewidth', 1.5, 'color', cpal(8,:))
        ylim([0, inf+1])
    end
end

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'isothermal_cstr/simulation/report_3', '-dpdf', '-r300')
