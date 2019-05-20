%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamical Modelling and Intelligent Control of Chemical Reactive Systems
% 
% Author: Otacilio Bezerra Leite Neto
%
%   report_codes.m
%       This script contains the codes used to generate the images displayed in the 
%       Chapter 4 for the monograph on report/ (Diagrams and Illustrations are not included)
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

%% Figure 4.2 %%
% - Simulation Parameters
% Linear Model Index
idx = 25;
% Time
t = (0:0.05:15.95)'; T = numel(t);
% Initial Conditions
ue = iso_cstr.oper.U(idx,:); xe = iso_cstr.oper.X(idx,:);

% Reference Signal
r = [ones(1,T/4)*xe(2) ones(1,T/4)*xe(2)-2 ones(1,T/4)*xe(2)+1 ones(1,T/4)*xe(2)-1];

% Disturbance signal
w = randn(T, 2) .* 0.03;      % Process Noise
z = randn(T, 1) .* 0.1;      % Measurement Noise

% Linear Model
A = iso_cstr.ss_model.A(idx);   B = iso_cstr.ss_model.B(idx);
C = iso_cstr.ss_model.C;        D = iso_cstr.ss_model.D;

iso_cstr.ss_model.C = [0 1];
iso_cstr.ss_model.D = [0];
iso_cstr.sizeY      =  1 ;

% Controller and Observer
Q = diag([30, 30]);
R = 10;

L = [0; 0];

% - Simulation of the Outputs
[~, yout, xout, uout] = simulate(iso_cstr, idx, t, r, xe, 'switch-lqr', Q, R, 100, L, w, z);

% - Visualization of the Simulation
figure(1);
subplot(1,2,1)
plot(t, uout, 'linewidth', 1.5, 'color', cpal(1,:)); 
xlabel("Time (min)"), ylabel("u = \Delta u + u_o")

subplot(2,2,2)
plot(t, xout(1,:), 'linewidth', 1, 'linestyle', '--', 'color', cpal(2,:)); hold on;
plot(t, yout(1,:), 'linewidth', 1.5, 'linestyle', '-', 'color', cpal(2,:)); hold on;
ylabel("Outflow Concentration - C_a (mol/l)")

subplot(2,2,4), 
plot(t, r, 'linestyle', '--', 'color', 'black'); hold on;
plot(t, xout(2,:), 'linewidth', 1, 'linestyle', '--', 'color', cpal(2,:)); hold on;
plot(t, yout(2,:), 'linewidth', 1.5, 'linestyle', '-', 'color', cpal(2,:)); hold on; 
xlabel("Time (min)"), ylabel("Outflow Concentration - C_b (mol/l)")
