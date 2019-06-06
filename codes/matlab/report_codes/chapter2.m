%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamical Modelling and Intelligent Control of Chemical Reactive Systems
% 
% Author: Otacilio Bezerra Leite Neto
%
%   report_ch2_codes.m
%       This script contains the codes used to generate the images displayed in the 
%       Chapter2 2for the monograph on report/ (Diagrams and Illustrations are not included)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preamble %%
% cd /home/minho/Documents/Minho/Control_ChemPlants/
% cd /home/minhotmog/Dropbox/Research/TCC/Codes/
clc; clear all; close all;
addpath("utils/")

% Sets the Default Renderer to tbe the Painters
set(0, 'DefaultFigureRenderer', 'painters');

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

%% Figure 2.3 %%
% - Simulation of the Nonlinear plat to obtain steady-state values
[~, y] = simulate(iso_cstr.sysVar, (0:.1:4)', ones(1,41)*3.03, [0 0 0 0]);
x_o = y(:,end);

% - Real simulation parameters
t = (0:0.01:5.99)';                                          
X_0 = x_o;
U = [ones(1, numel(t)/4)*1  ones(1, numel(t)/4)*3  ones(1, numel(t)/4)*2 ones(1, numel(t)/4)*4];

% - State-space model selection
A = iso_cstr.ss_model.A_full(x_o, 3.03);     B = iso_cstr.ss_model.B_full(x_o, 3.03);
C = eye(4);                                                  D = zeros(4,1);
iso_cstr_lin = ss(A, B, C, D);

% - Simulation of the Outputs
[~, y] = simulate(iso_cstr.sysVar, t, U+3.03, X_0);
[y_lin, ~, ~] = lsim(iso_cstr_lin, U, t, X_0-x_o);

% - Visualization of the Simulation
figure(1)
clf
[ha, pos] = tight_subplot(1,1,[.02 .1],[.1 .01],[.1 .01]);
axes(ha(1))
plot(t, U, 'k-', 'linewidth', 1)
ylim([0, 5]), ytickformat("%.1f")
xlabel("Time (min)"), ylabel("u (1/min)")

% - Exporting the Visualization to an Image
figname = "report_codes/figs/report_ch2_1_1";
fig = gcf; fig.PaperPositionMode = 'auto'; 
print('-bestfit', figname, '-dpdf', '-r300')
system("pdfcrop " + figname + ".pdf " + figname + ".pdf");

figure(2)
clf
[ha, pos] = tight_subplot(2,2,[.02 .1],[.1 .01],[.1 .01]);
for i = 1:4
    axes(ha(i))
    plot(t, y(i,:), 'linewidth', 1, 'color', cpal(5+i, :)); hold on;
    plot(t, y_lin(:,i)+x_o(i), 'linewidth', 1, 'linestyle', '--', 'color', cpal(5+i, :))
    ylabel(strcat("x_", int2str(i), " (mol/l)"))
end
axes(ha(1)); set(ha(1), "XTickLabel", []), ytickformat("%.1f")
axes(ha(2)); set(ha(2), "XTickLabel", [])
axes(ha(3)); xlabel("Time (min)")
axes(ha(4)); xlabel("Time (min)")

% - Exporting the Visualization to an Image
figname = "report_codes/figs/report_ch2_1_2";
fig = gcf; fig.PaperPositionMode = 'auto'; 
print('-bestfit', figname, '-dpdf', '-r300')
system("pdfcrop " + figname + ".pdf " + figname + ".pdf");

%% Figure 2.4 %%
% - Simulation of the Nonlinear plat to obtain steady-state values
[~, y] = simulate(iso_cstr.sysVar, (0:.1:4)', ones(1,41)*3.03, [0 0 0 0]);
x_o = y(1:2,end);

% - Real simulation parameters
t = (0:0.1:1)';      
U = ones(1,numel(t)) * 3;
U_zero = zeros(1,numel(t));

% - State-space model selection
A = iso_cstr.ss_model.A_l(x_o, 3.03);     B = iso_cstr.ss_model.B_l(x_o, 3.03);
C = eye(2);                               D = zeros(2,1);
iso_cstr_lin = ss(A, B, C, D);

% - Simulation of the Outputs
[y_lin, ~, ~]      = lsim(iso_cstr_lin, U, t, [10 10]/3.5);
[y_lin_nat, ~, ~]  = lsim(iso_cstr_lin, U_zero, t, [10 10]/3.5);
y_lin_forc = y_lin - y_lin_nat;

% - Visualization of the Simulation
figure(1); clf
[ha, pos] = tight_subplot(2,1,[.02 .1],[.1 .01],[.1 .01]);
axes(ha(1))
plot(t, y_lin(:,1), 'linewidth', 1, 'color', cpal(6, :)); hold on;
plot(t, y_lin_nat(:,1), 'linewidth', 1, 'linestyle', '--', 'color', cpal(6, :)); hold on;
plot(t, y_lin_forc(:,1), 'linewidth', 1, 'linestyle', ':', 'color', cpal(6, :)); hold on;
ylabel("\Delta x_1 (mol/l)")
set(ha(1), "XTickLabel", [])

axes(ha(2))
plot(t, zeros(1,numel(t)), 'k-'); hold on;
plot(t, y_lin(:,2), 'linewidth', 1, 'color', cpal(8, :)); hold on;
plot(t, y_lin_nat(:,2), 'linewidth', 1, 'linestyle', '--', 'color', cpal(8, :)); hold on;
plot(t, y_lin_forc(:,2), 'linewidth', 1, 'linestyle', ':', 'color', cpal(8, :)); hold on;
ylabel("\Delta x_2 (mol/l)")
xlabel("Time (min)")
ylim([-1, 3.5])

% - Exporting the Visualization to an Image
figname = "report_codes/figs/report_ch2_2_1";
fig = gcf; fig.PaperPositionMode = 'auto'; fig.PaperSize = [4 3];
print('-bestfit', figname, '-dpdf', '-r300')
system("pdfcrop " + figname + ".pdf " + figname + ".pdf");


figure(2); clf
[ha, pos] = tight_subplot(1,1,[.02 .1],[.1 .01],[.1 .01]);
axes(ha(1))
scatter(y_lin(end,1), y_lin(end,2), 'ko')
z = zeros(size(t));
surface([y_lin(:,1) y_lin(:,1)]',[y_lin(:,2) y_lin(:,2)]',[z z]',[t t]',...
        'facecol','no', 'edgecol','interp', 'linew',1.5); hold on
surface([y_lin_nat(:,1) y_lin_nat(:,1)]',[y_lin_nat(:,2) y_lin_nat(:,2)]',[z z]',[t t]',...
        'facecol','no', 'edgecol','interp', 'linew',1, 'LineStyle', '--'); hold on
surface([y_lin_forc(:,1) y_lin_forc(:,1)]',[y_lin_forc(:,2) y_lin_forc(:,2)]',[z z]',[t t]',...
        'facecol','no', 'edgecol','interp', 'linew',1, 'LineStyle', ':');
colormap jet
c = colorbar;
c.Label.String = 'Time (min)';

ylabel("\Delta x_2 (mol/l)")
xlabel("\Delta x_1 (mol/l)")
grid

% - Exporting the Visualization to an Image
figname = "report_codes/figs/report_ch2_2_2";
fig = gcf; fig.PaperPositionMode = 'auto'; fig.PaperSize = [4 4];
print('-bestfit', figname, '-dpdf', '-r300')
system("pdfcrop " + figname + ".pdf " + figname + ".pdf");


%% Figure 2.5 %%
% - Simulation of the Nonlinear plat to obtain steady-state values
[~, y] = simulate(iso_cstr.sysVar, (0:.1:4)', ones(1,41)*3.03, [0 0 0 0]);
x_o = y(1:2,end);

% - Real simulation parameters
t_s = (0:0.01:1.99)';      
U = ones(1,numel(t_s));
U_zero = zeros(1,numel(t_s));

% - State-space model selection
A = iso_cstr.ss_model.A_l(x_o, 3.03);

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
print('-bestfit', 'isothermal_cstr/simulation/report_ch2_3', '-dpdf', '-r300')

%% Figure 2.6 %%
% - Real simulation parameters
t_s = (0:0.1:49.9)';      
U = ones(1,numel(t_s));
U_zero = zeros(1,numel(t_s));

% - State-space model selection
A = [-0.1 0.5; -0.5 -0.1];

% - Calculation of the Matrix Exponential
syms t
e_At = expm(A * t);

% - Visualization of the Matrix Exponential
figure(4);
for i = 1:2
    for j = 1:2
        subplot(2,2,(i-1)*2+j)
        plot(t_s, exp(A(1,1).*t_s), 'linewidth', 1, 'linestyle', '--', 'color', [0.6 0.6 0.6]); hold on
        plot(t_s, -exp(A(1,1).*t_s), 'linewidth', 1, 'linestyle', '--', 'color', [0.6 0.6 0.6]); hold on
        plot(t_s, double(subs(e_At(i,j), t_s)), 'linewidth', 1.5, 'color', cpal(8,:))
        ylim([-inf, inf+1])
    end
end

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'isothermal_cstr/simulation/report_ch2_4', '-dpdf', '-r300')

%% Figure 2.5 %%
% - Simulation of the Nonlinear plat to obtain steady-state values
[~, y] = simulate(iso_cstr.sysVar, (0:.1:4)', ones(1,41)*3.03, [0 0 0 0]);
x_o = y(1:2,end);

% - Real simulation parameters
t1 = (0:0.01:1.99)';
t2 = (0:0.01:59.99)';
U1 = ones(1,numel(t1));
U2 = ones(1,numel(t2));

% - State-space model selection
A1 = iso_cstr.ss_model.A_l(x_o, 3.03);     B1 = iso_cstr.ss_model.B_l(x_o, 3.03);
C = [1 0];                                 D = 0;
iso_cstr_lin1 = ss(A1, B1, C, D);

A2 = [-0.1 0.5; -0.5 -0.1];     B2 = [1; 1];
iso_cstr_lin2 = ss(A2, B2, C, D);

% - Simulation of the Outputs
[y_lin1, ~, ~] = step(iso_cstr_lin1, t1);
[y_lin2, ~, ~] = step(iso_cstr_lin2, t2);

% - Visualization of the Simulation
figure(5);

subplot(1,2,1), plot(t1, U1 * y_lin1(end), 'linestyle', '--', 'linewidth', 1, 'color', [0.6 0.6 0.6]); hold on
plot(t1, y_lin1, 'linewidth', 1.5, 'color', cpal(8,:))
xlabel("Time"), ylabel("x^{(1)}_1(t)")

subplot(1,2,2), plot(t2, U2 * y_lin2(end), 'linestyle', '--', 'linewidth', 1, 'color', [0.6 0.6 0.6]); hold on
plot(t2, y_lin2, 'linewidth', 1.5, 'color', cpal(8,:))
xlabel("Time"), ylabel("x^{(2)}_1(t)")

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'isothermal_cstr/simulation/report_ch2_5', '-dpdf', '-r300')

%% Figure 2.6 %%
% - Real simulation parameters
t = (0:0.01:19.99)';
U = ones(1,numel(t));

% - State-space model selection
A1 = [-0.4];                    B1 = [0.137];
A2 = [-0.2];                      B2 = [0.069];
A3 = [-0.25 0.7; -0.7 -0.25];   B3 = [0.2; 0.2];
A4 = [-0.5   1;   -1 -0.5];     B4 = [0.287; 0.287];

iso_cstr_lin1 = ss(A1, B1, 1, 0);
iso_cstr_lin2 = ss(A2, B2, 1, 0);
iso_cstr_lin3 = ss(A3, B3, [1 0], 0);
iso_cstr_lin4 = ss(A4, B4, [1 0], 0);

% - Simulation of the Outputs
[y_lin1, ~, ~] = step(iso_cstr_lin1, t);
[y_lin2, ~, ~] = step(iso_cstr_lin2, t);
[y_lin3, ~, ~] = step(iso_cstr_lin3, t);
[y_lin4, ~, ~] = step(iso_cstr_lin4, t);

% - Visualization of the Simulation
figure(6);

subplot(1,2,1)
plot(real(eig(A1)), imag(eig(A1)), 'x', 'color', cpal(6,:), 'markersize', 12); hold on
plot(real(eig(A2)), imag(eig(A2)), 'x', 'color', cpal(7,:), 'markersize', 12); hold on
plot(real(eig(A3)), imag(eig(A3)), 'x', 'color', cpal(8,:), 'markersize', 12); hold on
plot(real(eig(A4)), imag(eig(A4)), 'x', 'color', cpal(9,:), 'markersize', 12)
xlabel("Real Axis"); ylabel("Imaginary Axis")
ylim([-1.5, 1.5]); xlim([-0.75, 0])
sgrid()

subplot(1,2,2),
plot(t, y_lin1, 'linewidth', 1.5, 'color', cpal(6,:)); hold on
plot(t, y_lin2, 'linewidth', 1.5, 'color', cpal(7,:)); hold on
plot(t, y_lin3, 'linewidth', 1.5, 'color', cpal(8,:)); hold on
plot(t, y_lin4, 'linewidth', 1.5, 'color', cpal(9,:))
xlabel("Time"), ylabel("Forced Response")
ylim([0, 0.6]);

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'isothermal_cstr/simulation/report_ch2_6', '-dpdf', '-r300')

%% Figure 2.7 %%

% - Defining the basis
P1 = [1 0; 0 1];
P2 = [3 2; 1 2]; slope1 = P2(1,1) / P2(2,1); slope2 = P2(1,2) / P2(2,2);

x1 = P1 * [1; 3];
x2 = pinv(P2) * [1; 3];

% - Visualization of the Simulation
figure(7);

subplot(1,2,1)
plot([5, -5], [0, 0], 'color', [0.5 0.5 0.5]); hold on
plot([0, 0], [-5, 5], 'color', [0.5 0.5 0.5]); hold on

quiver(0, 0, P1(1,1)+0.1, P1(2,1), 'linewidth', 1.5, 'color', [0.2 0.2 0.2]); hold on
quiver(0, 0, P1(1,2), P1(2,2)+0.1, 'linewidth', 1.5, 'color', [0.2 0.2 0.2]); hold on
quiver(0, 0, x1(1,1)+0.1, x1(2,1)+0.3, 'linewidth', 1.5, 'color', cpal(8,:));
xlim([-3.5, 3.5]); ylim([-1, 3.5])
grid()
set(gca, 'YTickLabel', [])
set(gca, 'XTickLabel', [])

subplot(1,2,2)
plot([5, -5], [0, 0], 'color', [0.5 0.5 0.5]); hold on 
plot([0, 0], [-5, 5], 'color', [0.5 0.5 0.5]); hold on

quiver(0, 0, P2(1,1)+0.3, P2(2,1)+0.1, 'linewidth', 1.5, 'color', [0.2 0.2 0.2]); hold on
quiver(0, 0, P2(1,2)+0.2, P2(2,2)+0.2, 'linewidth', 1.5, 'color', [0.2 0.2 0.2]); hold on
quiver(0, 0, x1(1,1)+0.1, x1(2,1)+0.3, 'linestyle', ':', 'linewidth', 1.5, 'color', cpal(8,:));
quiver(0, 0, x2(1,1)-0.1, x2(2,1)+0.2, 'linewidth', 1.5, 'color', cpal(8,:));
xlim([-3.5, 3.5]); ylim([-1, 3.5])
grid()
set(gca, 'YTickLabel', [])
set(gca, 'XTickLabel', [])

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'isothermal_cstr/simulation/report_ch2_7', '-dpdf', '-r300')

%% Figure 2.8 %%
% - Real simulation parameters
t = (0:0.01:24.99)';
U = ones(1,numel(t));

% - State-space model selection
A1 = [-0.1 1.5; -1.5 -0.1];   B1 = [1; 1];
A2 = [0 1.5; -1.5 0];          B2 = [1; 1];
A3 = [0.1 1.5; -1.5 0.1];     B3 = [1; 1];

iso_cstr_lin1 = ss(A1, B1, [1 0], 0);
iso_cstr_lin2 = ss(A2, B2, [1 0], 0);
iso_cstr_lin3 = ss(A3, B3, [1 0], 0);

% - Simulation of the Outputs
[y_lin1, ~, ~] = step(iso_cstr_lin1, t);
[y_lin2, ~, ~] = step(iso_cstr_lin2, t);
[y_lin3, ~, ~] = step(iso_cstr_lin3, t);

% - Visualization of the Simulation
figure(8);

subplot(1,2,1)
plot([5, -5], [0, 0], 'color', [0.5 0.5 0.5]); hold on
plot([0, 0], [-5, 5], 'color', [0.5 0.5 0.5]); hold on
plot(real(eig(A1)), imag(eig(A1)), 'x', 'color', cpal(6,:), 'markersize', 12); hold on
plot(real(eig(A2)), imag(eig(A2)), 'x', 'color', cpal(7,:), 'markersize', 12); hold on
plot(real(eig(A3)), imag(eig(A3)), 'x', 'color', cpal(8,:), 'markersize', 12); hold on
xlabel("Real Axis"); ylabel("Imaginary Axis")
ylim([-2.5, 2.5]); xlim([-0.2, 0.2])
grid()

subplot(1,2,2),
plot(t, y_lin1, 'linewidth', 1.5, 'color', cpal(6,:)); hold on
plot(t, y_lin2, 'linewidth', 1.5, 'color', cpal(7,:)); hold on
plot(t, y_lin3, 'linewidth', 1.5, 'color', cpal(8,:)); hold on
xlabel("Time"), ylabel("Forced Response")
ylim([-12 12])

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'isothermal_cstr/simulation/report_ch2_8', '-dpdf', '-r300')

%% Figure 2.9 %%
% - Real simulation parameters
t = (0:0.01:24.99)';
U = ones(1,numel(t));

% - Sinusoidal signals
y1 = 1 * cos(0.75*t + 0);
y2 = 1.2 * cos(0.75*t - 1);

% - Visualization of the Simulation
figure;
plot(t, y1, 'linewidth', 1.5, 'color', cpal(8,:))
title('u(t) = M_i cos(\omega t + \phi_i)')
ylim([-1.5, 1.5])
set(gca, 'YTickLabel', [])
set(gca, 'XTickLabel', [])

fig = gca;
print('isothermal_cstr/simulation/report_ch2_9_1', '-dpng', '-r300')

figure;
plot(t, y1, 'linewidth', 1.5, 'color', cpal(8,:), 'linestyle', '--'); hold on
plot(t, y2, 'linewidth', 1.5, 'color', cpal(8,:))
ylim([-1.5, 1.5])
set(gca, 'YTickLabel', [])
set(gca, 'XTickLabel', [])

fig = gca;
print('isothermal_cstr/simulation/report_ch2_9_2', '-dpng', '-r300')

%% Figure 2.10 %%
% - Simulation of the Nonlinear plat to obtain steady-state values
[~, y] = simulate(iso_cstr.sysVar, (0:.1:4)', ones(1,41)*3.03, [0 0 0 0]);
x_o = y(1:2,end);

% - State-space model selection
A = iso_cstr.ss_model.A_l(x_o, 3.03);   B = iso_cstr.ss_model.B_l(x_o, 3.03);

% - Visualization of the Simulation
figure(10);
subplot(1,2,1)
bodeplot(tf(ss(A, B, [1 0], 0)), 'b-'); hold on 
lineHandle = findobj(gcf,'Type','line','-and','Color','b');
set(lineHandle,'Color',cpal(8,:));

bodeplot(tf(ss(A, B, [0 1], 0)), 'b-');
lineHandle = findobj(gcf,'Type','line','-and','Color','b');
set(lineHandle,'Color',cpal(4,:));
grid()

subplot(2,2,2)
nyquistplot(tf(ss(A, B, [1 0], 0)), 'b-'), xlabel(""), xlim([-0.1, 0.75])
lineHandle = findobj(gcf,'Type','line','-and','Color','b');
set(lineHandle,'Color',cpal(8,:));

subplot(2,2,4)
nyquistplot(tf(ss(A, B, [0 1], 0)), 'b-'), title(""), xlim([-0.15, 0.025])
lineHandle = findobj(gcf,'Type','line','-and','Color','b');
set(lineHandle,'Color',cpal(4,:));

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'isothermal_cstr/simulation/report_ch2_10', '-dpdf', '-r300')