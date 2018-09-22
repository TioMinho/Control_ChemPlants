%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamical Modelling and Intelligent Control of Chemical Reactive Systems
% 
% Author: Otacilio Bezerra Leite Neto
%
%   This script contains the instructions for running and visualize experiments 
%   over the Chemical Reactive systems models.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preamble %%
cd /home/minhotmog/Dropbox/Research/TCC/Codes/
clc; clear all; close all;

% Sets the Default Rendere to tbe the Painters
set(0, 'DefaultFigureRenderer', 'painters');

% Some colors
cpal = [209 17 65;    % 1 - Metro Red
             0 177 89;     % 2 - Metro Green
            0 174 219;    % 3 - Metro Blue
           243 119 53;   % 4 - Metro Orange
           255 196 37;   % 5 - Metro Yellow
                              ]/255;  

%% %%%%%%%%%%%%
%  ISOTHERMAL CSTR $$
%  %%%%%%%%%%%%%%
%% Model Loading %%
run isothermal_cstr/iso_model.m

%% Model Non-Linear Simulation %%
% - Simulation Parameters
% Time
t = (0:0.1:15)';                                         
% Initial Conditions
X_0 = [0 0 0 0];                                     
% Input Signal
U = [ones(15,1)*0; ones(40,1); ones(numel(t)-55,1)*0];
 
% - Simulation of the Outputs
[~, y] = simulate(iso_cstr.model, t, U, X_0);
    
% - Visualization of the Simulation
figure(1);
subplot(1,2,1), plot(t, U, 'linewidth', 1.5), title("Input Signal")
xlabel("Time (min)"), ylabel("Input Flow-rate (m^3/min)")
grid()

subplot(1,2,2), p = plot(t, y); title("Output Signals")
p(1).LineWidth = 1.5; p(2).LineWidth = 1.5;
p(3).LineStyle='--'; p(4).LineStyle='--';
legend("C_A", "C_B", "C_C", "C_D") 
xlabel("Time (min)"), ylabel("Outflow Concentration (mol/l)")
grid()

% - Exporting the Visualization to an Image
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-bestfit', 'simulation/isoCSTR_sim_01', '-dpdf', '-r300')

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
figure;
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
print('-bestfit', 'simulation/exoCSTR_sim_01', '-dpdf', '-r300')