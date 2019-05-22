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

xe = reshape(exo_cstr.oper.X(13,13,:), [1 4]); ue = exo_cstr.oper.U(13, :);
oper.xe = xe';
oper.ue = ue';

exo_cstr.ss_model.C = [0 1 0 0; 0 0 1 0];
exo_cstr.ss_model.D = [0 0; 0 0];
exo_cstr.sizeY      =  2 ;

%% Figure 4.5 %%
% - Simulation Parameters
% Time
t = (0:0.1:0.5)'; T = numel(t);

% Reference Signal
r = [ones(2,T).*[xe(2); xe(3)], ...
     ones(2,T).*[xe(2); xe(3)].*[1.1;1.01], ...
     ones(2,T).*[xe(2); xe(3)].*[0.9;0.99], ...
     ones(2,T).*[xe(2); xe(3)]];

% Disturbance signal
w = randn(T, 4) .* [.01 .01 .01 .01];      % Process Noise
z = randn(T, 2) .* [.01 .01];      % Measurement Noise

% Linear Model
A = exo_cstr.ss_model.A_l(xe, ue);   B = exo_cstr.ss_model.B_l(xe, ue);
C = exo_cstr.ss_model.C;             D = exo_cstr.ss_model.D;

% Controller and Observer
Q = diag([25, 25, 1, 1, 10, 10]);
R = diag([10 25]);

% - Simulation of the Outputs
try
    [~, yout, xout, uout] = simulate(exo_cstr, oper, t, r, xe, 'lqgi', Q, R, numel(t), w, z);

    % Exports the Visualization to a PDF
    exo_sim.t = t; exo_sim.oper = oper; exo_sim.r = r; exo_sim.xe = xe;
    exo_sim.type = 'lqri'; exo_sim.Q = Q; exo_sim.R = R; 
    exo_sim.xout = xout; exo_sim.yout = yout; exo_sim.uout = uout;
   
    timeNow = datetime('now', 'TimeZone', 'local', 'Format', 'dMMMy_HHmmssZ');
    save("data/simulations/exoSim_"+char(timeNow), 'exo_sim')

    xout = xout + xe';
    yout = yout + xe';
    uout = uout + ue';

    % - Visualization of the Simulation
    figure(1);
    subplot(2,2,1)
    plot(t, uout(1,:), 'linewidth', 1.5, 'color', cpal(8,:));
    subplot(2,2,3)
    plot(t, uout(2,:), 'linewidth', 1.5, 'color', cpal(8,:));

    subplot(2,2,2) 
    plot(t, r(1,:), 'linestyle', '--', 'color', 'black'); hold on;
    plot(t, yout(2,:), 'linewidth', 1.5, 'color', cpal(8,:)); hold on;     
    xlabel("Time (min)")

    subplot(2,2,4) 
    plot(t, r(2,:), 'linestyle', '--', 'color', 'black'); hold on;
    plot(t, yout(3,:), 'linewidth', 1.5, 'color', cpal(8,:)); hold on;     
    xlabel("Time (min)")

    subplot(2,2,1), ylabel("u_1 = \Delta u_1 + u_{1o}")
    subplot(2,2,3), ylabel("u_2 = \Delta u_2 + u_{2o}")
    subplot(2,2,2), ylabel("x_1 = \Delta x_1 + x_{o1} (mol/L)")
    subplot(2,2,4), ylabel("x_2 = \Delta x_2 + x_{o2} (^o C)")

catch err
    fprintf("Erro!\n");
    rethrow(err)
    
end

%% !!!!!!!!!! AUTOMATIC EXPERIMENT MODE !!!!!!!!!!
run report_experiments/experiment_parameters.m

for i=1:n_experiments
    % - Simulation of the Outputs
    try
        [~, yout, xout, uout] = simulate(exo_cstr, oper, exp_param{i}.t, exp_param{i}.r, ...
                                         exp_param{i}.x_0, exp_param{i}.type, exp_param{i}.Q, ... 
                                         exp_param{i}.R, numel(exp_param{i}.t), exp_param{i}.w, exp_param{i}.z);

        % Exports the Visualization to a PDF
        xout = xout + xe';
        yout = yout + xe';
        uout = uout + ue';
        
        exp_param{i}.xout = xout; exp_param{i}.yout = yout; exp_param{i}.uout = uout;

        timeNow = datetime('now', 'TimeZone', 'local', 'Format', 'dMMMy_HHmmssZ');
        save("data/simulations/exoSim_"+char(timeNow), 'exp_param')
        
    catch err
        fprintf("Erro!\n");
        rethrow(err)

    end
end