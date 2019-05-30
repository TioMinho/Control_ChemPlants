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
        lines(30)*255
       ]/255;  

%% %%%%%%%%%%%%%%%
%  EXOTHERMAL CSTR  %
%%%%%%%%%%%%%%%%%%%%%
%% Model Loading %%
run exothermal_cstr/exo_model.m

xe = [1.235 0.9 134.14 128.95]; ue = [18.83 -4495.7];
oper.xe = xe';
oper.ue = ue';

exo_cstr.ss_model.C = [0 1 0 0; 0 0 1 0];
exo_cstr.ss_model.D = [0 0; 0 0];
exo_cstr.sizeY      =  2 ;

%% %%%%%%%%%%%%%%%
%  EXOTHERMAL CSTR  %
%%%%%%%%%%%%%%%%%%%%%
%%
% = Simulation Parameters = %
% - Time and Initial States
t = (0:0.005:2.695)'; T = numel(t);
x_0 = xe;

% - Reference Signal
r = [ones(1,50)*xe(2) xe(2) + 0.6*square(linspace(0, 3*pi-1, T-50)) + 0.4*square(linspace(pi, 4*pi-1, T-50)); xe(3)*ones(1,T); ];
%r = ones(2,T).*[xe(2); xe(3)];

% - Noise and Disturbance signals
w = mvnrnd([0; 0; 0; 0], diag([0.1 0.1 0.1 0.1]), T);  % Process Noise
z = mvnrnd([0; 0], diag([0.0001, 0.1]), T);            % Measurement Noise

W = [ 0*ones(1,T); 
      zeros(1,T)]'; %mvnrnd([0; 0; 0], diag([0.01, 0.01, 0.1]), T);

W(100:110,1) = 0.4;
W(300:310,1) = -0.4;

W(200:210,2) = -2;
W(400:410,2) = +2;
  
% = Controller and Estimator Definitions = %
controller.type = "lqri";

controller.oper.xe = xe'; controller.oper.ue = ue';
controller.N = T;

controller.Q = diag([1, 1, 1, 1, 1e3, 1e3]);
controller.R = diag([1, 1e-3]);

controller.Q_k = cov(w);
controller.R_k = cov(z);

% ----
Qi_values = [1e4 1e5 1e4 1e5];
Q_values = [10:5:500];
R_values = [10:5:500];


%% ---------------------------
% = Simulation of the Outputs = %
try
    for i=1:1
%     controller.Q = diag(Q_values(randperm(numel(Q_values), 4)));
%     controller.R = diag(R_values(randperm(numel(R_values), 2)));
    
%     controller.Q = diag([Q_values(randperm(numel(Q_values), 4)), Qi_values(1) 1e-4]); controller.Q
%     controller.R = diag(R_values(randperm(numel(R_values), 2))); controller.R
    
    [~, yout, xout, uout] = simulate(exo_cstr, t, r, x_0, "control", controller, "disturbance", W);

    % Exports the Visualization to a PDF
    exp_param.model = exo_cstr;
    exp_param.t = t; exp_param.r = r; exp_param.x_0 = x_0;
    exp_param.w = w; exp_param.z = z;
    exp_param.controller = controller;
    
    exp_param.xout = xout; 
    exp_param.yout = yout + xe'; 
    exp_param.uout = uout + ue';
   
    timeNow = datetime('now', 'TimeZone', 'local', 'Format', 'dMMMy_HHmmssZ');
    save("data/simulations/"+controller.type+"/exoSim_"+char(timeNow), 'exp_param')

    % - Visualization of the Simulation
    figure(1);
    subplot(3,2,1)
    plot(exp_param.t, exp_param.uout(1,:), 'linewidth', 1.25, 'color', cpal(10+i,:)); hold on
    subplot(3,2,2)
    plot(exp_param.t, exp_param.uout(2,:), 'linewidth', 1.25, 'color', cpal(10+i,:)); hold on
    xlabel("Time (hr)")
    
    subplot(3,2,3) 
    plot(exp_param.t, exp_param.r(1,:), 'linestyle', '--', 'color', [0.3 0.3 0.3]); hold on;
    %plot(exp_param.t, exp_param.yout(2,:), 'linestyle', '--', 'linewidth', 1.5, 'color', cpal(10+i,:)); hold on;
    plot(exp_param.t, exp_param.xout(2,:), 'linewidth', 1.25, 'color', cpal(10+i,:)); hold on;
%     s = scatter(exp_param.t, exp_param.xout(2,:)+exp_param.z(:,1)', 'x', 'MarkerEdgeColor', cpal(8,:)); hold on;     
%     s.MarkerFaceAlpha = 0.3;
%     s.MarkerEdgeAlpha = 0.3;
    ylabel("Concentration - \rho_B (^o C)")

    subplot(3,2,5) 
    %plot(exp_param.t, exp_param.yout(1,:), 'linestyle', '--', 'linewidth', 1.5, 'color', cpal(10+i,:)); hold on;
    plot(exp_param.t, exp_param.xout(1,:), 'linewidth', 1.25, 'color', cpal(10+i,:)); hold on;
    ylabel("Concentration - \rho_A (^o C)")
    xlabel("Time (hr)")
    
    subplot(3,2,4) 
    plot(exp_param.t, exp_param.r(2,:), 'linestyle', '--', 'color', [0.3 0.3 0.3]); hold on;
    %plot(exp_param.t, exp_param.yout(3,:), 'linestyle', '--', 'linewidth', 1.5, 'color', cpal(10+i,:)); hold on;
    plot(exp_param.t, exp_param.xout(3,:), 'linewidth', 1.25, 'color', cpal(10+i,:)); hold on;     
%     s = scatter(exp_param.t, exp_param.xout(3,:)+exp_param.z(:,2)', 'x', 'MarkerEdgeColor', cpal(8,:)); hold on;     
%     s.MarkerFaceAlpha = 0.3;
%     s.MarkerEdgeAlpha = 0.3;
    ylabel("Temperature - Reactor (^o C)")

    subplot(3,2,6) 
    %plot(exp_param.t, exp_param.yout(4,:), 'linestyle', '--', 'linewidth', 1.5, 'color', cpal(10+i,:)); hold on;
    plot(exp_param.t, exp_param.xout(4,:), 'linewidth', 1.25, 'color', cpal(10+i,:)); hold on;     
    ylabel("Temperature - Coolant (^o C)")
    xlabel("Time (hr)")

    subplot(3,2,1)
    plot(exp_param.t, ones(1,T)*35, 'r:', 'linewidth', 1.25), hold on
    plot(exp_param.t, ones(1,T)*5, 'r:', 'linewidth', 1.25)
    ylabel("Input Flow-Rate (1/hr)")
    ylim([4, 36])
    
    subplot(3,2,2)
    plot(exp_param.t, ones(1,T)*0, 'r:', 'linewidth', 1.25), hold on
    plot(exp_param.t, ones(1,T)*-8500, 'r:', 'linewidth', 1.25)
    ylabel("Cooling Capacity (kJ/hr)")
    ylim([min(exp_param.uout(2,:)), max(exp_param.uout(2,:))])
    %ylim([-8700, 200])
    
    subplot(3,2,1), hold off
    subplot(3,2,2), hold off
    subplot(3,2,3), hold off
    subplot(3,2,4), hold off
    subplot(3,2,5), hold off
    subplot(3,2,6), hold off
    
    end
    
catch err
    fprintf("Erro!\n");
    rethrow(err)
    
end

timeNow = datetime('now', 'TimeZone', 'local', 'Format', 'dMMMy_HHmmssZ');
figname = "report_experiments/figs/exoSim_control_" + char(timeNow);
fig = gcf; fig.PaperPositionMode = 'auto'; 
print('-bestfit', figname, '-dpdf', '-r300')
system("pdfcrop " + figname + ".pdf " + figname + ".pdf");


%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% !!!!!!!!!! AUTOMATIC EXPERIMENT MODE !!!!!!!!!!
while(true)
    tic
    % - Simulation of the Outputs
    try
        exp_param = generateRandomExperiment(xe);
        [~, yout, xout, uout] = simulate(exo_cstr, oper, exp_param.t, exp_param.r, ...
                                         exp_param.x_0, exp_param.type, exp_param.Q, ... 
                                         exp_param.R, numel(exp_param.t), exp_param.w, exp_param.z);

        % Exports the Visualization to a PDF
        exp_param.xout = xout + xe';
        exp_param.yout = yout + xe';
        exp_param.uout = uout + ue';
        
        timeNow = datetime('now', 'TimeZone', 'local', 'Format', 'dMMMy_HHmmssZ');
        simulation = exp_param;
        save("data/simulations/exoSim_"+char(timeNow), 'simulation')
        
    catch err
        fprintf("======= Erro! ======\n");
        disp(getReport(err,'extended'));

    end
    toc
end

%%
% - Visualization of the Simulation
foldername = "data/simulations/lqri/";
folder = dir(foldername);
files = {folder.name}; files = files(130:end);
for name = files
    try
    clear exp_param exo_sim simulation
    load(foldername+name)
    
    if(exist("exp_param", "var"))
        simulation = exp_param;
        save(foldername+name, 'simulation')
    elseif(exist("exo_sim", "var"))
        simulation = exo_sim;
        save(foldername+name, 'simulation')
    end

    simulation
    diag(simulation.Q)
    
    figure(2);
    subplot(2,2,1)
    plot(simulation.t, simulation.uout(1,:), 'linewidth', 1.5, 'color', cpal(8,:));
    subplot(2,2,3)
    plot(simulation.t, simulation.uout(2,:), 'linewidth', 1.5, 'color', cpal(8,:));

    subplot(2,2,2) 
    if(size(simulation.r,1) > 0), plot(simulation.t, simulation.r(1,:), 'linestyle', '--', 'color', 'black'); hold on;
    end
    if(strcmp(simulation.type, 'lqg') || strcmp(simulation.type, 'lqgi'))
        s = scatter(simulation.t, simulation.xout(2,:), 'x', 'MarkerEdgeColor', cpal(8,:)); hold on;     
        s.MarkerFaceAlpha(0.5)
    end
    plot(simulation.t, simulation.yout(2,:), 'linewidth', 1.5, 'color', cpal(8,:)); hold on;     
    xlabel("Time (min)")

    subplot(2,2,4) 
    if(size(simulation.r,1) > 1), plot(simulation.t, simulation.r(2,:), 'linestyle', '--', 'color', 'black'); hold on;
    end
    if(strcmp(simulation.type, 'lqg') || strcmp(simulation.type, 'lqgi'))
        s = scatter(simulation.t, simulation.xout(3,:), 'x', 'MarkerEdgeColor', cpal(8,:)); hold on;     
        s.MarkerFaceAlpha(0.5)
    end
    plot(simulation.t, simulation.yout(3,:), 'linewidth', 1.5, 'color', cpal(8,:)); hold on;     
    xlabel("Time (min)")

    subplot(2,2,1), ylabel("u_1 = \Delta u_1 + u_{1o}"), hold off
    subplot(2,2,3), ylabel("u_2 = \Delta u_2 + u_{2o}"), hold off
    subplot(2,2,2), ylabel("x_1 = \Delta x_1 + x_{o1} (mol/L)"), hold off
    subplot(2,2,4), ylabel("x_2 = \Delta x_2 + x_{o2} (^o C)"), hold off
    
    pause
    end
end


%%
foldername = "data/simulations/";
folder = dir(foldername);
files = {folder.name}; files = files(3:end);
for name = files
    clear exp_param exo_sim simulation
    load(foldername+name)
    save(foldername+exp_param.type+"/"+name, 'exp_param')
end