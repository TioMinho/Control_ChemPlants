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

%% Figure 4.5 %%
% - Simulation Parameters
% Linear Model Index
idx = 25;
% Time
t = (0:0.05:59.95)'; T = numel(t);
% Initial Conditions
ue = iso_cstr.oper.U(idx,:); xe = iso_cstr.oper.X(idx,:);
controller.oper.ue = ue; controller.oper.xe = xe';
x_0 = xe;

% Reference Signal
r = [ones(1,T/4)*xe(2) ones(1,T/4)*xe(2)-0.2 ones(1,T/4)*xe(2)+0.2 ones(1,T/4)*xe(2)];

% Noise signal
w = mvnrnd([0; 0], diag([.001 .001]), T);
z = mvnrnd([0], diag([0.001]), T);      

% Linear Model
A = iso_cstr.ss_model.A(idx);   B = iso_cstr.ss_model.B(idx);
C = iso_cstr.ss_model.C;        D = iso_cstr.ss_model.D;

iso_cstr.ss_model.C = [0 1];
iso_cstr.ss_model.D = [0];
iso_cstr.sizeY      =  1 ;

% Controller and Observer
controller.type = "lqgi";
controller.Q = diag([10, 10, 1]);
controller.R = 1e-1; 

controller.Q_k = diag(diag(cov(w)));
controller.R_k = diag(diag(cov(z))); 

% - Simulation of the Outputs
[~, yout, xout, uout] = simulate(iso_cstr, t, r, x_0, "control", controller, 'noise', [w z]);
xout = xout + xe'; 
yout = yout + xe'; 
uout = uout + ue';

% - Visualization of the Simulation
figure(1);clf
[ha, pos] = tight_subplot(2,2,[.03 .1],[.1 .01],[.1 .01]);

axes(ha(1))
plot(t, uout(1,:), 'linewidth', 1, 'color', [0.2 0.2 0.2]); hold on
xlabel("Time (min)")
ylabel("u = \Delta u + u_o (1/min)")

axes(ha(2))
plot(t, r(1,:), 'linewidth', 1, 'linestyle', '--', 'color', [0.3 0.3 0.3]); hold on;
plot(t, yout(2,:), 'linestyle', '--', 'linewidth', 1, 'color', cpal(6,:)); hold on;
s = scatter(t, xout(2,:)+z(:,1)', '.', 'MarkerEdgeColor', cpal(6,:)); hold on;     
s.MarkerFaceAlpha = 0.4;
s.MarkerEdgeAlpha = 0.4;
plot(t, xout(2,:), 'linewidth', 1, 'color', cpal(6,:)); hold on;
set(ha(2), "XTickLabel", [])
ylabel("x_2 = \Delta x_2 + x_{o2} (mol/l)")

axes(ha(4))
plot(t, yout(1,:), 'linestyle', '--', 'linewidth', 1, 'color', cpal(6,:)); hold on;
plot(t, xout(1,:), 'linewidth', 1, 'color', cpal(6,:)); hold on;
xlabel("Time (min)")
ylabel("x_1 = \Delta x_1 + x_{o1} (mol/l)")

axes(ha(3)), set(gca,'Visible','off')

% - Exporting the Visualization to an Image
figname = "report_codes/figs/report_ch4_1";
fig = gcf; fig.PaperPositionMode = 'auto'; fig.PaperSize = [10 5];
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