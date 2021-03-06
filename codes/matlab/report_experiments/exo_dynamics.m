%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamical Modelling and Intelligent Control of Chemical Reactive Systems
% 
% Author: Otacilio Bezerra Leite Neto
%
%   exo_dynamics.m
%       This script contains the instructions for running and visualize experiments 
%       over the dynamics of the Exothermal Continuous-Stirred Tank (CSTR) system.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preamble %%
clc; close all; clear all;
% cd /home/minho/Documents/Minho/Control_ChemPlants/
% cd /home/minhotmog/Dropbox/Research/TCC/Codes/

% Sets the Default Rendere to tbe the Painters
set(0, 'DefaultFigureRenderer', 'painters');
addpath("utils")

% Some colors
cpal = [209 17 65;    % 1 - Metro Red
        0 177 89;     % 2 - Metro Green
        0 174 219;    % 3 - Metro Blue
        243 119 53;   % 4 - Metro Orange
        255 196 37;   % 5 - Metro Yellow
        %%
        217,83,79;     % 6 - Bootstrap Red
        91,192,222;    % 7 - Bootstrap Light Blue
        92,184,92;     % 8 - Bootstrap Green
        66,139,202;    % 9 - Bootstrap Blue
        255,167,0;     % 10 - Google Yellow
        lines(10)*255
       ]/255;  
                          
%% %%%%%%%%%%%%
%  EXOTHERMAL CSTR %%
%  %%%%%%%%%%%%%%
%% Model Loading %%
run models/exo_model.m
%load('data/exo_cstr_model.mat')

xe = [1.235 0.9 134.14 128.95]; ue = [18.83 -4495.7];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%
%  Model Non-Linear Simulation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Parameters
% Time
t = (0:0.01:1)'; T = numel(t); sT = floor(T/3);
% Initial Conditions
X_0 = [xe 2.4192 0.2730];
% Input Signal
U = [ones(1,sT)*ue(1) ones(1,sT)*ue(1)*1.6 ones(1,T-2*sT)*ue(1)
     ones(1,sT)*ue(2) ones(1,sT)*ue(2)*1.6 ones(1,T-2*sT)*ue(2)];   
 
% Simulation of the Outputs
[~, y] = simulate(exo_cstr.sysVar, t, U, X_0);
    
% Visualization
figure(1); clf; fs = 18;
[ha, pos] = tight_subplot(2,2,[.03 .1],[.1 .01],[.1 .01]);

axes(ha(1)), plot(t, U(1,:), 'linewidth', 1, 'color', [0.3 0.3 0.3])
ylabel("Flow-rate (1/hr)")
ylim([min(U(1,:))-0.5, max(U(1,:))+0.5])
set(ha(1), 'XTickLabel', [])
% set(ha(1), "FontSize", fs)

axes(ha(2))
plot(t, y(1, :), 'linewidth', 1, 'color', cpal(11,:)); hold on
plot(t, y(2, :), 'linewidth', 1, 'color', cpal(12,:)); hold on
plot(t, y(5, :), 'linestyle', '--', 'linewidth', 1, 'color', cpal(13,:)); hold on
plot(t, y(6, :), 'linestyle', '--', 'linewidth', 1, 'color', cpal(14,:));
ylabel("Concentration (mol/l)")
set(ha(2),'XTickLabel', [])
% set(ha(2), "FontSize", fs)

axes(ha(3)), plot(t, U(2,:), 'linewidth', 1, 'color', [0.3 0.3 0.3])
xlabel("Time (hr)"), ylabel("Cooling Capacity (kJ/hr)")
ylim([min(U(2,:))-500, max(U(2,:))+500])
% set(ha(3), "FontSize", fs)

axes(ha(4)), 
plot(t, y(3, :), 'linewidth', 1, 'color', cpal(17,:)); hold on
plot(t, y(4, :), 'linewidth', 1, 'color', cpal(16,:));
xlabel("Time (hr)"), ylabel("Temperatures (�C)")
% set(ha(4), "FontSize", fs)

% - Exporting the Visualization to an Image
timeNow = datetime('now', 'TimeZone', 'local', 'Format', 'dMMMy_HHmmssZ');
figname = "report_experiments/figs/exoSim_dynamics_" + char(timeNow);
fig = gcf; 
fig.PaperUnits = "centimeters"; 
xSize = 1600; ySize = 600; xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
fig.Position = [xLeft yTop xSize ySize];
print('-bestfit', figname, '-dpdf', '-r300')
system("pdfcrop " + figname + ".pdf " + figname + ".pdf");

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Model Non-Linear Simulation - Disturbed %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Parameters
% Time
t = (0:0.01:1)'; T = numel(t);
% Initial Conditions
X_0 = [squeeze(exo_cstr.oper.X(10,10,:))' 0 0];
% Input Signal
U = [ones(1,T)*exo_cstr.oper.U(10,1)
     ones(1,T)*exo_cstr.oper.U(10,2)];   
 
W = [sin(t)'*0.1; sin(t)'*0.6; sin(t)'];
 
% Simulation of the Outputs
[~, y] = simulate(exo_cstr.sysVar_W, t, [U; W], X_0);
    
% Visualization
figure(1); clf
[ha, pos] = tight_subplot(2,2,[.03 .1],[.1 .01],[.1 .01]);

axes(ha(1)), plot(t, U(1,:), 'linewidth', 1.5, 'color', [0.3 0.3 0.3])
ylabel("Flow-rate (1/hr)")
set(ha(1), "XTickLabel", []) 

axes(ha(2))
p = plot(t, y(1, :), 'linewidth', 1.5, 'color', cpal(11,:)); hold on
p = plot(t, y(2, :), 'linewidth', 1.5, 'color', cpal(12,:)); hold on
p = plot(t, y(5, :), 'linestyle', '--', 'linewidth', 1, 'color', cpal(13,:)); hold on
p = plot(t, y(6, :), 'linestyle', '--', 'linewidth', 1, 'color', cpal(14,:));
xlabel("Time (hr)"), ylabel("Concentration (mol/l)")

axes(ha(3)), plot(t, U(2,:), 'linewidth', 1.5, 'color', [0.3 0.3 0.3])
ylabel("Cooling Capacity (kJ/hr)")
ylim([-7000 -4000])
set(ha(3), "XTickLabel", []) 

axes(ha(4))
plot(t, y(3, :), 'linewidth', 1.5, 'color', cpal(17,:)); hold on
plot(t, y(4, :), 'linewidth', 1.5, 'color', cpal(16,:));
xlabel("Time (hr)"), ylabel("Temperatures (ºC)")

% - Exporting the Visualization to an Image
timeNow = datetime('now', 'TimeZone', 'local', 'Format', 'dMMMy_HHmmssZ');
figname = "report_experiments/figs/exoSim_dynamics_" + char(timeNow);
fig = gcf; fig.PaperPositionMode = 'auto'; 
print('-bestfit', figname, '-dpdf', '-r300')
system("pdfcrop " + figname + ".pdf " + figname + ".pdf");


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Model Linearized Simulation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Simulation Parameters
% Time
t = (0:0.01:3.99)'; T = numel(t);

% Initial Conditions
U_ss = [18.83 -4495.7]; X_ss = [1.235 0.9 134.14 128.95];

% Input Signal
U = [ones(1,T-6*floor(T/7))*U_ss(1) ones(1,floor(T/7))*U_ss(1)*0.4 ones(1,floor(T/7))*U_ss(1) ones(1,floor(T/7))*U_ss(1) ones(1,floor(T/7))*U_ss(1) ones(1,floor(T/7))*U_ss(1)*0.4 ones(1,floor(T/7))*U_ss(1);
     ones(1,T-6*floor(T/7))*U_ss(2) ones(1,floor(T/7))*U_ss(2) ones(1,floor(T/7))*U_ss(2) ones(1,floor(T/7))*U_ss(2)*0.4 ones(1,floor(T/7))*U_ss(2) ones(1,floor(T/7))*U_ss(2)*0.4 ones(1,floor(T/7))*U_ss(2)];

U = U + randn(2,T) .* (0.0*U_ss');
 
% Linear Model
A = exo_cstr.ss_model.A_l(X_ss, U_ss); B = exo_cstr.ss_model.B_l(X_ss, U_ss);
C = exo_cstr.ss_model.C;               D = exo_cstr.ss_model.D;

exo_cstr_lin = ss(A, B, C, D);

% - Simulation of the Outputs
[~, y] = simulate(exo_cstr.sysVar, t, U, [X_ss 2.4192 0.2730]);
[y_lin, ~, ~] = lsim(exo_cstr_lin, U-U_ss', t);

y_lin = y_lin+X_ss;

%% Visualization
figure(1); clf
[ha, pos] = tight_subplot(2,2,[.03 .1],[.1 .01],[.1 .01]);

axes(ha(1)), plot(t, U(1,:), 'linewidth', 1, 'color', [0.3 0.3 0.3])
ylabel("Flow-rate (1/hr)")
ylim([min(U(1,:))-1, max(U(1,:))+1])
set(ha(1), "XTickLabel", [])

axes(ha(2))
plot(t, y(1, :), 'linestyle', '--', 'linewidth', 0.5, 'color', cpal(11,:)); hold on
plot(t, y(2, :), 'linestyle', '--', 'linewidth', 0.5, 'color', cpal(12,:)); hold on
plot(t, y_lin(:, 1), 'linestyle', '-', 'linewidth', 1, 'color', cpal(11,:)); hold on
plot(t, y_lin(:, 2), 'linestyle', '-', 'linewidth', 1, 'color', cpal(12,:));
ylabel("Concentration (mol/l)")
set(ha(2), "XTickLabel", [])

axes(ha(3)), plot(t, U(2,:), 'linewidth', 1, 'color', [0.3 0.3 0.3])
ylabel("Cooling Capacity (kJ/hr)")
ylim([min(U(2,:))-250, max(U(2,:))+250])
xlabel("Time (hr)")

axes(ha(4)), 
plot(t, y(3, :), 'linestyle', '--', 'linewidth', 0.5, 'color', cpal(17,:)); hold on
plot(t, y(4, :), 'linestyle', '--', 'linewidth', 0.5, 'color', cpal(16,:)); hold on
plot(t, y_lin(:, 3), 'linestyle', '-', 'linewidth', 1, 'color', cpal(17,:)); hold on
plot(t, y_lin(:, 4), 'linestyle', '-', 'linewidth', 1, 'color', cpal(16,:));
xlabel("Time (hr)"), ylabel("Temperatures (^oC)")

% - Exporting the Visualization to an Image
timeNow = datetime('now', 'TimeZone', 'local', 'Format', 'dMMMy_HHmmssZ');
figname = "report_experiments/figs/exoSim_dynamics_" + char(timeNow);
fig = gcf; 
fig.PaperUnits = "centimeters"; 
xSize = 1600; ySize = 600; xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
fig.Position = [xLeft yTop xSize ySize];
print('-bestfit', figname, '-dpdf', '-r300')
system("pdfcrop " + figname + ".pdf " + figname + ".pdf");

%% Visualization - Animated
timeNow = datetime('now', 'TimeZone', 'local', 'Format', 'dMMMy_HHmmssZ');
figname = "report_experiments/figs/exoSim_dynamics_" + char(timeNow);
vidfile = VideoWriter(figname + ".avi", "Motion JPEG AVI");

figure(1); clf

[ha, pos] = tight_subplot(2,2,[.05 .1],[.1 .01],[.1 .01]);
z = zeros(size(t));
open(vidfile);
for l = 1:numel(t)
    l_inf = max(1, l-numel(t)*0.15);

    axes(ha(1)), plot(t(1:l), U(1,1:l), 'linewidth', 1, 'color', [0.3 0.3 0.3])
    ylabel("Flow-rate (1/hr)")
    ylim([min(U(1,:))-1, max(U(1,:))+1])
    xlabel("Time (hr)")
    xlim([t(1), t(end)])

    axes(ha(2)), plot(t(1:l), U(2,1:l), 'linewidth', 1, 'color', [0.3 0.3 0.3])
    ylabel("Cooling Capacity (kJ/hr)")
    ylim([min(U(2,:))-250, max(U(2,:))+250])
    xlabel("Time (hr)")
    xlim([t(1), t(end)])
    
    axes(ha(3))
    scatter(xe(1), xe(2), 'ko'), hold on
    surface([y(1,l_inf:l)' y(1,l_inf:l)']',[y(2,l_inf:l)' y(2,l_inf:l)']',[z(l_inf:l) z(l_inf:l)]',[t(l_inf:l) t(l_inf:l)]',...
             'edgecol','interp', 'linew', 1); hold on
    surface([y_lin(l_inf:l,1) y_lin(l_inf:l,1)]',[y_lin(l_inf:l,2) y_lin(l_inf:l,2)]',[z(l_inf:l) z(l_inf:l)]',[t(l_inf:l) t(l_inf:l)]',...
             'edgecol','interp', 'linew',1, 'LineStyle', ':'); 
    colormap jet
    c = colorbar("NorthOutside");
    c.Label.String = 'Time (hr)';

    ylabel("\Delta x_2 (mol/l)")
    xlabel("\Delta x_1 (mol/l)")
    
    xlim([min(y(1,:))-0.1 max(y(1,:))+0.1])
    ylim([min(y(2,:))-0.05 max(y(2,:))+0.05])
    
    
    axes(ha(4))
    scatter(xe(3), xe(4), 'ko'), hold on
    surface([y(3,l_inf:l)' y(3,l_inf:l)']',[y(4,l_inf:l)' y(4,l_inf:l)']',[z(l_inf:l) z(l_inf:l)]',[t(l_inf:l) t(l_inf:l)]',...
             'edgecol','interp', 'linew', 1); hold on
    surface([y_lin(l_inf:l,3) y_lin(l_inf:l,3)]',[y_lin(l_inf:l,4) y_lin(l_inf:l,4)]',[z(l_inf:l) z(l_inf:l)]',[t(l_inf:l) t(l_inf:l)]',...
             'edgecol','interp', 'linew',1, 'LineStyle', '--'); 
    colormap jet
    c = colorbar("NorthOutside");
    c.Label.String = 'Time (hr)';

    ylabel("\Delta x_2 (mol/l)")
    xlabel("\Delta x_1 (mol/l)")
    
    xlim([min(y(3,:))-1 max(y(3,:))+1])
    ylim([min(y(4,:))-1 max(y(4,:))+1])
    
    %
    drawnow
    F(l) = getframe(gcf); 
    writeVideo(vidfile,F(l));
     
    axes(ha(1)), cla
    axes(ha(2)), cla
    axes(ha(3)), cla
    axes(ha(4)), cla
end

close(vidfile);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Steady-State Points Visualization %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3); clf
surface(exo_cstr.oper.U(:,1), exo_cstr.oper.U(:,2), exo_cstr.oper.X(:,:,2)', 'EdgeColor','none'); hold on
scatter3(ue(1), ue(2), 2, 200, [0.9 1 0.9], 'x', 'linewidth', 2)
ylim([-8500, 0]), ylabel("Cooling Capacity (kJ/h)"), xlabel("Flow-rate (1/hr)")
colorbar('NorthOutside')
colormap jet

% - Exporting the Visualization to an Image
timeNow = datetime('now', 'TimeZone', 'local', 'Format', 'dMMMy_HHmmssZ');
figname = "report_experiments/figs/exoSim_dynamics_" + char(timeNow);
fig = gcf; fig.PaperPositionMode = 'auto'; 
print('-bestfit', figname, '-dpdf', '-r300')
system("pdfcrop " + figname + ".pdf " + figname + ".pdf");

figure(4); clf
surface(exo_cstr.oper.U(:,1), exo_cstr.oper.U(:,2), exo_cstr.oper.X(:,:,3)', 'EdgeColor','none'); hold on
scatter3(ue(1), ue(2), 155, 200, [0.9 1 0.9], 'x', 'linewidth', 2)
ylim([-8500, 0]), ylabel("Cooling Capacity (kJ/h)"), xlabel("Flow-rate (1/hr)")
grid('off')
colorbar('NorthOutside')
colormap jet

% - Exporting the Visualization to an Image
timeNow = datetime('now', 'TimeZone', 'local', 'Format', 'dMMMy_HHmmssZ');
figname = "report_experiments/figs/exoSim_dynamics_" + char(timeNow);
fig = gcf; fig.PaperPositionMode = 'auto'; 
print('-bestfit', figname, '-dpdf', '-r300')
system("pdfcrop " + figname + ".pdf " + figname + ".pdf");

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Time Response Analysis (Matrix Exponential) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the system modes in a symbolic expression
syms t
[V, J] = jordan(A);
e_At = V*expm(J*t)*pinv(V);

% Simulate the contribution of each mode
t = 0:0.001:1.5;

%% Visualization
figure(5);

for j = 1:4 
for k = 1:4
    tic
    subplot(4,4,j+4*(k-1)), plot(t, subs(e_At(j,k), t), 'linewidth', 1.5, 'color', cpal(11,:))
    toc
end
end

% - Exporting the Visualization to an Image
timeNow = datetime('now', 'TimeZone', 'local', 'Format', 'dMMMy_HHmmssZ');
figname = "report_experiments/figs/exoSim_dynamics_" + char(timeNow);
fig = gcf; fig.PaperPositionMode = 'auto'; 
print('-bestfit', figname, '-dpdf', '-r300')
system("pdfcrop " + figname + ".pdf " + figname + ".pdf");

%% %%%%%%%%%%%%%%%%%%%%%%%%%%
%  Linearized Model Properties %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stability2
stabMatrix = zeros(50,50);
for i=1:50
    for j=1:50
        A = exo_cstr.ss_model.A_l(exo_cstr.oper.X(i,j,:), [exo_cstr.oper.U(i,1) exo_cstr.oper.U(j,2)]);
        stabMatrix(i,j) = max(eig(A));
    end
end

figure(6);
surface(exo_cstr.oper.U(:,1), exo_cstr.oper.U(:,2), stabMatrix');
ylim([-8500, 0]), ylabel("Cooling Capacity (kJ/h)"), xlabel("Flow-rate (1/hr)")
grid('off')
colorbar('NorthOutside')
colormap bone
shading interp

% - Exporting the Visualization to an Image
timeNow = datetime('now', 'TimeZone', 'local', 'Format', 'dMMMy_HHmmssZ');
figname = "report_experiments/figs/exoSim_dynamics_" + char(timeNow);
fig = gcf; fig.PaperPositionMode = 'auto'; 
print('-bestfit', figname, '-dpdf', '-r300')
system("pdfcrop " + figname + ".pdf " + figname + ".pdf");

%% Controllability
ctrbMatrix = zeros(50,50);
for i=1:50
    for j=1:50
        A = exo_cstr.ss_model.A_l(exo_cstr.oper.X(i,j,:), [exo_cstr.oper.U(i,1) exo_cstr.oper.U(j,2)]);
        B = exo_cstr.ss_model.B_l(exo_cstr.oper.X(i,j,:), [exo_cstr.oper.U(i,1) exo_cstr.oper.U(j,2)]);
        ctrbMatrix(i,j) = rank(ctrb(A,B));
    end
end

figure(7);
surface(exo_cstr.oper.U(:,1), exo_cstr.oper.U(:,2), ctrbMatrix');
ylim([-8500, 0]), ylabel("Cooling Capacity (kJ/h)"), xlabel("Flow-rate (1/hr)")
grid('off')
colorbar('NorthOutside')
colormap bone
shading interp

% - Exporting the Visualization to an Image
timeNow = datetime('now', 'TimeZone', 'local', 'Format', 'dMMMy_HHmmssZ');
figname = "report_experiments/figs/exoSim_dynamics_" + char(timeNow);
fig = gcf; fig.PaperPositionMode = 'auto'; 
print('-bestfit', figname, '-dpdf', '-r300')
system("pdfcrop " + figname + ".pdf " + figname + ".pdf");

%% Observability
obsvMatrix = zeros(50,50);
for i=1:50
    for j=1:50
        A = exo_cstr.ss_model.A_l(exo_cstr.oper.X(i,j,:), [exo_cstr.oper.U(i,1) exo_cstr.oper.U(j,2)]);
        C = [0 1 0 0; 0 0 1 0];
        obsvMatrix(i,j) = rank(obsv(A,C));
    end
end

figure(8);
surface(exo_cstr.oper.U(:,1), exo_cstr.oper.U(:,2), obsvMatrix');
ylim([-8500, 0]), ylabel("Cooling Capacity (kJ/h)"), xlabel("Flow-rate (1/hr)")
grid('off')
colorbar('NorthOutside')
colormap bone
shading interp

% - Exporting the Visualization to an Image
timeNow = datetime('now', 'TimeZone', 'local', 'Format', 'dMMMy_HHmmssZ');
figname = "report_experiments/figs/exoSim_dynamics_" + char(timeNow);
fig = gcf; fig.PaperPositionMode = 'auto'; 
print('-bestfit', figname, '-dpdf', '-r300')
system("pdfcrop " + figname + ".pdf " + figname + ".pdf");