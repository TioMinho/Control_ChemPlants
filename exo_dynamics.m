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
%  EXOTHERMAL CSTR $$
%  %%%%%%%%%%%%%%
%% Model Loading %%
run exothermal_cstr/exo_model.m

%% Model Non-Linear Simulation %%
% Simulation Parameters
% Time
t = (0:0.1:15)';
% Initial Conditions
X_0 = [0 0 85 79.8 0 0];
% Input Signal
U = [ones(numel(t),1)*0.55     ones(numel(t),1)*0];   
 
% Simulation of the Outputs
[~, y] = simulate(exo_cstr.model, t, U, X_0);
    
% Visualization
figure(1);
subplot(2,2,1), plot(t, U(:,1), 'linewidth', 1.5), title("Input Signal")
xlabel("Time (min)"), ylabel("Flow-rate (m^3/min)")
grid()

subplot(2,2,2), p = plot(t, y(:, [1 2 5 6])); title("Output Signals")
p(1).LineWidth = 1.5; p(2).LineWidth = 1.5;
p(3).LineStyle='--'; p(4).LineStyle='--';
legend("C_A", "C_B", "C_C", "C_D"),
xlabel("Time (min)"), ylabel("Concentration (mol/m^3)")
grid()

subplot(2,2,3), plot(t, U(:,2), 'linewidth', 1.5)
xlabel("Time (min)"), ylabel("Cooling Capacity (MJ/min)")
grid()

subplot(2,2,4), plot(t, y(:, 3:4), 'linewidth', 1.5);
legend("T", "T_C"), 
xlabel("Time (min)"), ylabel("Temperatures (K)")
grid()

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'exothermal_cstr/simulation/exoCSTR_sim_01', '-dpdf', '-r300')

%% Model Linearized Simulation %%
% - Simulation Parameters
% Linear Model Index
idx = [13 13];
% Time
t = (0:0.1:10)';                                         
% Initial Conditions
U_0 = [exo_cstr.oper.U(idx(1),1) exo_cstr.oper.U(idx(2),2)]; X_0 = squeeze(exo_cstr.oper.X(idx(1),idx(2),:));
% Input Signal
U = [sin(2*t)*0.1     sin(2*t)*-0.05];

% Linear Model
A = exo_cstr.ss_model.A( (idx(1)-1)*25 + idx(2) );
B = exo_cstr.ss_model.B( (idx(1)-1)*25 + idx(2) );
exo_cstr_lin = ss(A, B, exo_cstr.ss_model.C, exo_cstr.ss_model.D);

% - Simulation of the Outputs
[~, y] = simulate(exo_cstr.model, t, U_0+U, [X_0' 0 0]);
[y_lin, ~, ~] = lsim(exo_cstr_lin, U, t);

% - Visualization of the Simulation
figure(2);
subplot(2,3,1), plot(t, U_0(1)+U(:,1), 'linewidth', 1.5), title("Input Signal")
ylim([5/60, 35/60]), xlabel("Time (min)"), ylabel("Flow-rate (m^3/min)")
grid()

subplot(2,3,2), p = plot(t, y(:, [1 2 5 6])); title("Non-Linear CSTR")
p(1).LineWidth = 1.5; p(2).LineWidth = 1.5;
p(3).LineStyle='--'; p(4).LineStyle='--';
legend("C_A", "C_B", "C_C", "C_D"),
ylim([0 inf]), xlabel("Time (min)"), ylabel("Concentration (mol/m^3)")
grid()

subplot(2,3,3), plot(t, y_lin(:, 1:2)+X_0(1:2)', 'linewidth', 1.5); title("Linearized CSTR")
legend("C_A", "C_B"),
ylim([0 inf]), xlabel("Time (min)"), ylabel("Concentration (mol/m^3)")
grid()

subplot(2,3,4), plot(t, U_0(2)+U(:,2), 'linewidth', 1.5)
ylim([-8.5/60, 0]), xlabel("Time (min)"), ylabel("Cooling Capacity (MJ/min)")
grid()

subplot(2,3,5), plot(t, y(:, [3 4]), 'linewidth', 1.5) 
legend("T", "T_C"), 
xlabel("Time (min)"), ylabel("Concentration (mol/m^3)")
grid()

subplot(2,3,6), plot(t, y_lin(:, 3:4)+X_0(3:4)', 'linewidth', 1.5);
legend("T", "T_C"), 
xlabel("Time (min)"), ylabel("Temperatures (K)")
grid()

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'exothermal_cstr/simulation/exoCSTR_sim_lin_01', '-dpdf', '-r300')

%% Operation Points Visualization %%
figure(3);
subplot(1,2,1), surf(exo_cstr.oper.U(:,1), exo_cstr.oper.U(:,2), exo_cstr.oper.X(:,:,2)', ...
                                    ccmap.area(:,:,1:3), 'AlphaData', ccmap.area(:,:,4), 'FaceAlpha', 'flat')
title("Operation Points - C_B")
xlabel("Input Flow-rate (m^3/min)"), ylabel("Cooling Capacity (MJ/h)"), zlabel("Outflow Concentration (mol/m^3)")

subplot(1,2,2), surf(exo_cstr.oper.U(:,1), exo_cstr.oper.U(:,2), exo_cstr.oper.X(:,:,3)', ...
                                    ccmap.area(:,:,1:3), 'AlphaData', ccmap.area(:,:,4), 'FaceAlpha', 'flat')
title("Operation Points - T")
xlabel("Input Flow-rate (m^3/min)"), ylabel("Cooling Capacity (MJ/h)"), zlabel("Tank Temperature (K)")

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'exothermal_cstr/simulation/exoCSTR_sim_op_01', '-dpdf', '-r300')

%% System Modes Visualization %%
t = 0:0.1:4;
figure(4);
for i = 1:exo_cstr.oper.size
    tic
    pos = [mod(i,25)+(mod(i,25) == 0)*25 ceil(i/25)];
    
    subplot(2,2,1), plot(t, subs(exo_cstr.modes(i,1), t), 'color', ccmap.area(pos(1), pos(2), :)), 
    title("First Mode"), xlabel("Time (min)"), ylabel("e^{\lambda_1 t}"), hold on
    
    pseudo_mode = subs(exo_cstr.modes(i,2), t);
    subplot(2,2,2), plot(t, pseudo_mode, 'color', ccmap.area(pos(1), pos(2), :))
    title("Second Mode"), xlabel("Time (min)"), ylabel("e^{\alpha_2 t} cos(\omega_2 t)"), hold on
    
    subplot(2,2,3), plot(t, pseudo_mode, 'color', ccmap.area(pos(1), pos(2), :))
    title("Third Mode"), xlabel("Time (min)"), ylabel("e^{\alpha_3 t}cos(\omega_3 t)"), hold on
    
    subplot(2,2,4), plot(t, subs(exo_cstr.modes(i,4), t), 'color', ccmap.area(pos(1), pos(2), :))
    title("Fourth Mode"), xlabel("Time (min)"), ylabel("e^{\lambda_4 t}"), hold on
    
    drawnow
    fprintf("%d\n", i);
    toc
end

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'exothermal_cstr/simulation/exoCSTR_sim_modes_01', '-dpdf', '-r300')

%% Transfer Matrices Visualization %%
t = 0:0.1:4;
figure(5);
for i = 1:exo_cstr.oper.size
    tic
    pos = [mod(i,25)+(mod(i,25) == 0)*25 ceil(i/25)];
     
    for j = 1:4 
    for k = 1:4
        subplot(4,4,j+4*(k-1)), plot(t, subs(exo_cstr.trf_matrix(j,k,i), t), 'color', ccmap.area(pos(1), pos(2), :)), hold on
    end
    end
    
    drawnow
    fprintf("%d\n", i);
    toc
end

subplot(4,4,1), title("X_1"), ylabel("X_1")
subplot(4,4,2), title("X_2"), subplot(4,4,3), title("X_3"), subplot(4,4,4), title("X_4")
subplot(4,4,5), ylabel("X_2"), subplot(4,4,9), ylabel("X_3")
subplot(4,4,13), ylabel("X_4"), xlabel("Time (min)")
subplot(4,4,14), xlabel("Time (min)"), subplot(4,4,15), xlabel("Time (min)"), subplot(4,4,16), xlabel("Time (min)")

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'exothermal_cstr/simulation/exoCSTR_sim_expAt_01', '-dpdf', '-r300')

%% Stability Analysis - Pole-Zero Mapping %%
figure(6);
for i = 1:exo_cstr.oper.size
    tic
    pos = [mod(i,25)+(mod(i,25) == 0)*25 ceil(i/25)];
    
    pole = exo_cstr.poles(i,:);
    plot(real(pole), imag(pole), 'x', 'color', ccmap.area(pos(1), pos(2),:), 'markersize', 10); hold on
    drawnow
    toc
end
title("Pole-Zero Map")
xlabel("Real Axis (seconds^{-1})")
ylabel("Imaginary Axis (seconds^{-1})")
sgrid()

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'exothermal_cstr/simulation/exoCSTR_stab_pz_01', '-dpdf', '-r300')

%% Stability Analysis - Bode Plots %%
figure(7);
for i = 1:exo_cstr.oper.size
    tic
    pos = [mod(i,25)+(mod(i,25) == 0)*25 ceil(i/25)];
    
    bodeplot(ss(exo_cstr.ss_model.A(i), exo_cstr.ss_model.B(i), exo_cstr.ss_model.C, exo_cstr.ss_model.D), 'b'), hold on;

    lineHandle = findobj(gcf,'Type','line','-and','Color','b');
    set(lineHandle,'Color',ccmap.area(pos(1), pos(2),:));
    
    drawnow
    toc
end
grid()

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'exothermal_cstr/simulation/exoCSTR_stab_bode_01', '-dpdf', '-r300')

%% Stability Analysis - Nyquist Plots %%
figure(8);
for i = 1:iso_cstr.oper.size
    tic
    pos = [mod(i,25)+(mod(i,25) == 0)*25 ceil(i/25)];
    
    nyquistplot(ss(exo_cstr.ss_model.A(i), exo_cstr.ss_model.B(i), exo_cstr.ss_model.C, exo_cstr.ss_model.D), 'b'), hold on;

    lineHandle = findobj(gcf,'Type','line','-and','Color','b');
    set(lineHandle,'Color',ccmap.area(pos(1), pos(2),:));
    
    drawnow
    toc
end
grid()

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'exothermal_cstr/simulation/exoCSTR_stab_nyquist_01', '-dpdf', '-r300')

