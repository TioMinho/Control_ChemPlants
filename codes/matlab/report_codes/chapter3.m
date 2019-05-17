%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamical Modelling and Intelligent Control of Chemical Reactive Systems
% 
% Author: Otacilio Bezerra Leite Neto
%
%   report_codes.m
%       This script contains the codes used to generate the images displayed in the 
%       Chapter 3 for the monograph on report/ (Diagrams and Illustrations are not included)
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

%% Figure 3.2 %%
% - Simulation of the Nonlinear plat to obtain steady-state values
u_o = 3.03;
[~, y] = simulate(iso_cstr.sysVar, (0:.1:4)', ones(1,41)*u_o, [0 0 0 0]);
x_o = y(1:2,end);

% - Real simulation parameters
t = (0:0.01:2)';      
r = zeros(1,numel(t));

% - State-space model selection
A = iso_cstr.ss_model.A_l(x_o, 3.03);     B = iso_cstr.ss_model.B_l(x_o, 3.03);
C = eye(2);                               D = zeros(2,1);

% - Controllers
K = {place(A,B,[-4, -1.5]); place(A,B,[-8, -5]); place(A,B,[-4+3j, -4-3j])};

% - Closed-Loop Systems
A_cl = {ss(A-B*K{1}, zeros(size(B)), C, D), ss(A-B*K{2}, zeros(size(B)), C, D), ss(A-B*K{3}, zeros(size(B)), C, D)};

% - Perform the simulation
Delta_x = [lsim(A_cl{1}, r, t, -x_o), lsim(A_cl{2}, r, t, -x_o), lsim(A_cl{3}, r, t, -x_o)];

% - Visualization of the Simulation
figure(1);
subplot(2,2,2), plot(t, r+x_o(1), 'k--'), hold on
subplot(2,2,4), plot(t, r+x_o(2), 'k--'),  hold on

for i = 1:3
    subplot(1,2,1)
    plot(t, (r - K{i}*Delta_x(:, 2*i-1:2*i)')+ u_o, 'linewidth', 1.5, 'color', cpal(5+i, :)); hold on;
    
    subplot(2,2,2)
    plot(t, Delta_x(:,2*i-1)+x_o(1), 'linewidth', 1.5, 'color', cpal(5+i, :)); hold on;
    
    subplot(2,2,4)
    plot(t, Delta_x(:,2*i)+x_o(2), 'linewidth', 1.5, 'color', cpal(5+i, :)); hold on;
end

subplot(1,2,1), ylabel("u = (\Delta u + u_{o})"), xlabel("Time (s)"), legend(["K_1", "K_2", "K_3"])
subplot(2,2,2), ylabel("x_1 = (\Delta x_1 + x_{o1}) (mol/l)")
subplot(2,2,4), ylim([-inf, x_o(2)+1]), ylabel("x_2 = (\Delta x_2 + x_{o2}) (mol/l)"), xlabel("Time (s)")

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'report_codes/figs/report_ch3_1', '-dpdf', '-r300')