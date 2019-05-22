%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamical Modelling and Intelligent Control of Chemical Reactive Systems
% 
% Author: Otacilio Bezerra Leite Neto
%
%   exo_control.m
%       This script contains the instructions for running and visualize experiments 
%       over the control of the Exothermic Continuous-Stirred Tank (CSTR) system.
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
run exothermal_cstr/exo_model.m

xe = [1.235 0.9 134.14 128.95]'; ue = [18.83 -4495.7]';
oper.xe = xe;
oper.ue = ue;

%% Figure 4.5 %%
% - Simulation Parameters
% Linear Model Index
idx = 13;
% Time
t = (0:0.5:19.5)'; T = numel(t);

% Reference Signal
r = [ones(2,T/4).*[xe(2); xe(3)]  ones(2,T/4).*[xe(2); xe(3)].*[1.1 1]' ones(2,T/4).*[xe(2); xe(3)].*[0.9 1]' ones(2,T/4).*[xe(2); xe(3)]];

% Linear Model
A = exo_cstr.ss_model.A_l(xe, ue);   B = exo_cstr.ss_model.B_l(xe, ue);
C = exo_cstr.ss_model.C;             D = exo_cstr.ss_model.D;

exo_cstr.ss_model.C = [0 1 0 0; 0 0 1 0];
exo_cstr.ss_model.D = [0 0; 0 0];
exo_cstr.sizeY      =  2 ;

% Controller and Observer
Q = diag([5, 5, 5, 5, 1e12, 1e5]);
R = diag([1 1]);

% - Simulation of the Outputs
[~, yout, xout, uout] = simulate(exo_cstr, oper, t, r, xe', 'lqri', Q, R, numel(t));
xout = xout + xe;
yout = yout + xe;
uout = uout + ue;

% - Visualization of the Simulation
figure(1);
subplot(2,2,1)
plot(t, uout(1,:), 'linewidth', 1.5, 'color', cpal(8,:));
subplot(2,2,3)
plot(t, uout(2,:), 'linewidth', 1.5, 'color', cpal(8,:));

subplot(2,2,2) 
plot(t, r(1,:), 'linestyle', '--', 'color', 'black'); hold on;
plot(t, yout(2,:), 'linewidth', 1.5, 'marker', 'o', 'linestyle', 'none', 'color', cpal(8,:)); hold on;     
xlabel("Time (min)")

subplot(2,2,4) 
plot(t, r(2,:), 'linestyle', '--', 'color', 'black'); hold on;
plot(t, yout(3,:), 'linewidth', 1.5, 'marker', 'o', 'linestyle', 'none', 'color', cpal(8,:)); hold on;     
xlabel("Time (min)")

subplot(2,2,1), ylabel("u = \Delta u + u_o")
subplot(2,2,3), ylabel("x_2 = \Delta x_2 + x_{o2} (mol/l)")