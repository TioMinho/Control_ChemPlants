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
addpath("utils")

% Some colors
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
run models/iso_model.m
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
figure(1); clf
[ha, pos] = tight_subplot(2,2,[.02 .1],[.1 .01],[.1 .01]);

axes(ha(2)), plot(t, r+x_o(1), 'k--'), hold on
axes(ha(4)), plot(t, r+x_o(2), 'k--'),  hold on

for i = 1:3
    axes(ha(1))
    plot(t, (r - K{i}*Delta_x(:, 2*i-1:2*i)')+ u_o, 'linewidth', 1, 'color', cpal(5+i, :)); hold on;
    
    axes(ha(2))
    plot(t, Delta_x(:,2*i-1)+x_o(1), 'linewidth', 1, 'color', cpal(5+i, :)); hold on;
    
    axes(ha(4))
    plot(t, Delta_x(:,2*i)+x_o(2), 'linewidth', 1, 'color', cpal(5+i, :)); hold on;
end

axes(ha(1)), ylabel("u = (\Delta u + u_{o}) (1/min)"), xlabel("Time (min)"), legend(["K_1", "K_2", "K_3"])
axes(ha(2)), ylabel("x_1 = (\Delta x_1 + x_{o1}) (mol/l)"), ytickformat("%.1f"), set(gca, "XTickLabel", []);
axes(ha(4)), ylim([-inf, x_o(2)+1]), ylabel("x_2 = (\Delta x_2 + x_{o2}) (mol/l)"), xlabel("Time (min)")
axes(ha(3)), set(gca,'Visible','off')

% - Exporting the Visualization to an Image
figname = "report_codes/figs/report_ch3_1";
fig = gcf; fig.PaperPositionMode = 'auto'; 
print('-bestfit', figname, '-dpdf', '-r300')
system("pdfcrop " + figname + ".pdf " + figname + ".pdf");

%% Figure 3.4 %%
% - Simulation of the Nonlinear plat to obtain steady-state values
u_o = 3.03;
[~, y] = simulate(iso_cstr.sysVar, (0:.1:4)', ones(1,41)*u_o, [0 0 0 0]);
x_o = y(1:2,end);

% - State-space model selection
A = iso_cstr.ss_model.A_l(x_o, 3.03);     B = iso_cstr.ss_model.B_l(x_o, 3.03);
C = [0, 1];                               D = 0;

ny = size(C,1); nu = size(B,2); nx = size(A,1);

% - Real simulation parameters
t = (0:0.01:11.99)'; T = numel(t);
r = [ones(1,T/4)*x_o(2) ones(1,T/4)*x_o(2)-0.2 ones(1,T/4)*x_o(2)+0.2 ones(1,T/4)*x_o(2)-0.4];

% - Closed-Loop Tracking Systems
K = {[place(A,B,[-4, -1.5]), -25]; [place(A,B,[-8, -5]), -25]; [place(A,B,[-4+3j, -4-3j]), -50]};
A_cl = {[A+B*K{1}(1:nx) B*K{1}(nx+1:end); -C zeros(ny)], ...
        [A+B*K{2}(1:nx) B*K{2}(nx+1:end); -C zeros(ny)], ...
        [A+B*K{3}(1:nx) B*K{3}(nx+1:end); -C zeros(ny)]};
      
B_cl = [zeros(nx, ny); eye(ny)];
C_i = [C, zeros(ny, ny)];
D_i = zeros(ny, ny);

% - Closed-Loop Systems
ia_sys = {ss(A_cl{1}, B_cl, eye(nx+ny), D_i), ... 
          ss(A_cl{2}, B_cl, eye(nx+ny), D_i), ...
          ss(A_cl{3}, B_cl, eye(nx+ny), D_i)};

% - Perform the simulation
Delta_x = [lsim(ia_sys{1}, r, t, -[x_o*0; r(1)-x_o(2)]), ...
           lsim(ia_sys{2}, r, t, -[x_o*0; r(1)-x_o(2)]), ...
           lsim(ia_sys{3}, r, t, -[x_o*0; r(1)-x_o(2)])];

% - Visualization of the Simulation
figure(1); clf
[ha, pos] = tight_subplot(1,2,[.02 .1],[.1 .01],[.1 .01]);

axes(ha(2)), plot(t, r+x_o(1), 'k--'), hold on

for i = 1:3
    axes(ha(1))
    plot(t, (r - K{i}*Delta_x(:, 3*i-2:3*i)')+ u_o, 'linewidth', 1, 'color', cpal(5+i, :)); hold on;
    
    axes(ha(2))
    plot(t, Delta_x(:,3*i-1)+x_o(1), 'linewidth', 1, 'color', cpal(5+i, :)); hold on;
end

axes(ha(1)), ylabel("u = (\Delta u + u_{o})"), xlabel("Time (min)"), legend(["K_1", "K_2", "K_3"])
axes(ha(2)), ylabel("x_2 = (\Delta x_2 + x_{o2}) (mol/l)"), xlabel("Time (min)")

% - Exporting the Visualization to an Image
figname = "report_codes/figs/report_ch3_2";
fig = gcf; fig.PaperPositionMode = 'auto'; 
print('-bestfit', figname, '-dpdf', '-r300')
system("pdfcrop " + figname + ".pdf " + figname + ".pdf");

%% Figure 3.6 %%
% - Equivalent Input-Output System and Closed-Loop representation
G = tf([.5 1.5],[1 2.5  1.5 0]);
T = feedback(2*G,1);

% - Calculate the Bode Plot Margins
margin(T)
grid("off")

lineHandle = findobj(gcf,'Type','line','-and','Color',[0 0.4470 0.7410]);
set(lineHandle, "linewidth", 1.5);
title("")

% - Exporting the Visualization to an Image
figname = "report_codes\figs\report_ch3_3_1";
fig = gcf; fig.PaperPositionMode = 'auto'; 
print('-bestfit', figname, '-dpdf', '-r300')
system("pdfcrop " + figname + ".pdf " + figname + ".pdf");

% - Calculate the Bode Plot Margins
nyquist(T)
grid("off")

lineHandle = findobj(gcf,'Type','line','-and','Color',[0 0.4470 0.7410]);
set(lineHandle, "linewidth", 1.5);
title("")

% - Exporting the Visualization to an Image
figname = "report_codes\figs\report_ch3_3_2";
fig = gcf; fig.PaperPositionMode = 'auto'; 
print('-bestfit', figname, '-dpdf', '-r300')
system("pdfcrop " + figname + ".pdf " + figname + ".pdf");

