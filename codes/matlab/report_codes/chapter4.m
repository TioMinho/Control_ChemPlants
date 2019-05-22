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

%% Figure 4.5 %%
% - Simulation Parameters
% Linear Model Index
idx = 25;
% Time
t = (0:0.1:19.9)'; T = numel(t);
% Initial Conditions
ue = iso_cstr.oper.U(idx,:); xe = iso_cstr.oper.X(idx,:);

% Reference Signal
r = [ones(1,T/4)*xe(2) ones(1,T/4)*xe(2)-0.5 ones(1,T/4)*xe(2)+0.8 ones(1,T/4)*xe(2)+0.4];

% Disturbance signal
w = randn(T, 2) .* [0.5 0.05];      % Process Noise
z = randn(T, 1) .* 0.1;      % Measurement Noise

% Linear Model
A = iso_cstr.ss_model.A(idx);   B = iso_cstr.ss_model.B(idx);
C = iso_cstr.ss_model.C;        D = iso_cstr.ss_model.D;

iso_cstr.ss_model.C = [0 1];
iso_cstr.ss_model.D = [0];
iso_cstr.sizeY      =  1 ;

% Controller and Observer
Q = cell(2); R = cell(2);
Q{1} = diag([20, 20, 5e3]); Q{2} = diag([1, 1, 1e4]);
R{1} = 75; R{2} = 20;

L = [0; 0];

% - Simulation of the Outputs
xout = cell(2); yout = cell(2); uout = cell(2);
for i=1:2
    [~, yout{i}, xout{i}, uout{i}] = simulate(iso_cstr, idx, t, r, xe, 'lqgi', Q{i}, R{i}, numel(t), L, w, z);
    xout{i} = xout{i} + xe';
    yout{i} = yout{i} + xe';
    uout{i} = uout{i} + ue;
end

% - Visualization of the Simulation
figure(1);
for i=1:2
    subplot(2,2,i)
    plot(t, uout{i}, 'linewidth', 1.5, 'color', cpal(5+i,:));

    subplot(2,2,2+i) 
    plot(t, r, 'linestyle', '--', 'color', 'black'); hold on;
    s = scatter(t, xout{i}(2,:), 'x', 'MarkerEdgeColor', cpal(5+i,:)); s.MarkerEdgeAlpha = 0.6; hold on;
    plot(t, yout{i}(2,:), 'linewidth', 1.5, 'linestyle', '-', 'color', cpal(5+i,:)); hold on;     
    xlabel("Time (min)")
end

subplot(2,2,1), ylabel("u = \Delta u + u_o")
subplot(2,2,3), ylabel("x_2 = \Delta x_2 + x_{o2} (mol/l)")

% - Exporting the Visualization to an Image
figname = "report_codes/figs/report_ch4_1";
fig = gcf; fig.PaperPositionMode = 'auto'; 
print('-bestfit', figname, '-dpdf', '-r300')
system("pdfcrop " + figname + ".pdf " + figname + ".pdf");

%% Figure 4.5 %%
% - Simulation Parameters
% - Equivalent Input-Output System and Closed-Loop representation
G = tf(poly([+1 +2]), poly([-2 -1.25 3]));
T = feedback(G, 1);

figure(1)
nyquistplot(-9*G, 'k-')
lineHandle = findobj(gcf,'Type','line','-and', 'Color', [0 0 0]);
set(lineHandle, "linewidth", 1.5);
set(lineHandle, "color", cpal(8,:));

circle(-1, 0, 1, [.7 .7 .7])
ylim([-3 3])
xlim([-3 3])

title("")

% - Exporting the Visualization to an Image
figname = "report_codes/figs/report_ch4_2";
fig = gcf; fig.PaperPositionMode = 'auto'; 
print('-bestfit', figname, '-dpdf', '-r300')
system("pdfcrop " + figname + ".pdf " + figname + ".pdf");

% - Filled Circle Function
function h = circle(x,y,r,c)
    hold on
    th = 0:pi/100:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    fill(xunit, yunit, c)
    alpha(.4)
    hold off
end