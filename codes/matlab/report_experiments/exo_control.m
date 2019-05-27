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
        lines(10)*255
       ]/255;  

%% %%%%%%%%%%%%
%  EXOTHERMAL CSTR 
%  %%%%%%%%%%%%%%
%% Model Loading %%
run exothermal_cstr/exo_model.m

xe = [1.235 0.9 134.14 128.95]; ue = [18.83 -4495.7];
oper.xe = xe';
oper.ue = ue';

exo_cstr.ss_model.C = [0 1 0 0; 0 0 1 0];
exo_cstr.ss_model.D = [0 0; 0 0];
exo_cstr.sizeY      =  2 ;

%% Figure 4.5 %%
% - Simulation Parameters
% Time and Initial States
exp_param.t = (0:0.1:11.9)'; exp_param.T = numel(exp_param.t);
exp_param.x_0 = xe;
exp_param.oper = oper;

% Reference Signal
exp_param.r = [ones(2,floor(exp_param.T/6)).*[xe(2); xe(3)], ...
               ones(2,floor(exp_param.T/3)).*[xe(2); xe(3)].*[1.2;1], ...
               ones(2,floor(exp_param.T/3)).*[xe(2); xe(3)].*[0.8;1], ...
               ones(2,exp_param.T-5*floor(exp_param.T/6)).*[xe(2); xe(3)]];
%exp_param.r = ones(2,exp_param.T).*[xe(2); xe(3)];


% Disturbance signal
exp_param.w = randn(exp_param.T, 4) .* [.0 .0 .0 .0];      % Process Noise
exp_param.z = randn(exp_param.T, 2) .* [.0 .0];              % Measurement Noise

% Controller and Observer
exp_param.type = "lqri";

exp_param.Q = diag([20, 20, 25, 25]);
exp_param.R = diag([50 10]);

exp_param.N = exp_param.T;

Qi       = [12 11 10 9];
Q_values = [1 5:5:500];
R_values = [1 5:5:500];

w_values = [linspace(0.01, 0.1); linspace(0.01, 0.1); linspace(0.01, 0.03); linspace(0.01, 0.03)]';
z_values = [linspace(0.01, 0.1); linspace(0.01, 0.1)]';

% ---------------------------
% - Simulation of the Outputs
try
    for i=1:4
    exp_param.Q = diag([Q_values(randperm(numel(Q_values), 4)) 10^Qi(i) 10^Qi(i)]);
    exp_param.R = diag(R_values(randperm(numel(R_values), 2)));
    
    [~, yout, xout, uout] = simulate(exo_cstr, oper, exp_param.t, exp_param.r, exp_param.x_0, exp_param.type, ...
                                     exp_param.Q, exp_param.R, exp_param.N, exp_param.w, exp_param.z);

    % Exports the Visualization to a PDF
    exp_param.xout = xout; 
    exp_param.yout = yout + xe'; 
    exp_param.uout = uout + ue';
   
    timeNow = datetime('now', 'TimeZone', 'local', 'Format', 'dMMMy_HHmmssZ');
    save("data/simulations/exoSim_"+char(timeNow), 'exp_param')

    % - Visualization of the Simulation
    figure(1);
    subplot(2,2,1)
    plot(exp_param.t, exp_param.uout(1,:), 'linewidth', 1.5, 'color', cpal(9+i,:)); hold on
    subplot(2,2,3)
    plot(exp_param.t, exp_param.uout(2,:), 'linewidth', 1.5, 'color', cpal(9+i,:)); hold on
    xlabel("Time (hr)")
    
    subplot(2,2,2) 
    plot(exp_param.t, exp_param.r(1,:), 'linestyle', '--', 'color', [0.3 0.3 0.3]); hold on;
    plot(exp_param.t, exp_param.yout(2,:), 'linewidth', 1.5, 'color', cpal(9+i,:)); hold on;
%     s = scatter(exp_param.t, exp_param.xout(2,:), 'x', 'MarkerEdgeColor', cpal(8,:)); hold on;     
%     s.MarkerFaceAlpha = 0.3;
%     s.MarkerEdgeAlpha = 0.3;
    
    subplot(2,2,4) 
    plot(exp_param.t, exp_param.r(2,:), 'linestyle', '--', 'color', [0.3 0.3 0.3]); hold on;
    plot(exp_param.t, exp_param.yout(3,:), 'linewidth', 1.5, 'color', cpal(9+i,:)); hold on;     
%     s = scatter(exp_param.t, exp_param.xout(3,:), 'x', 'MarkerEdgeColor', cpal(8,:)); hold on;     
%     s.MarkerFaceAlpha = 0.3;
%     s.MarkerEdgeAlpha = 0.3;
    xlabel("Time (hr)")

    subplot(2,2,1)
    plot(exp_param.t, ones(1,exp_param.T)*35, 'r:', 'linewidth', 1.5), hold on
    plot(exp_param.t, ones(1,exp_param.T)*5, 'r:', 'linewidth', 1.5)
    ylabel("u_1 = \Delta u_1 + u_{1o}")
    ylim([4, 36])
    
    subplot(2,2,3)
    plot(exp_param.t, ones(1,exp_param.T)*0, 'r:', 'linewidth', 1.5), hold on
    plot(exp_param.t, ones(1,exp_param.T)*-8500, 'r:', 'linewidth', 1.5)
    ylabel("u_2 = \Delta u_2 + u_{2o}")
    ylim([-8700, 200])
    
    subplot(2,2,2), ylabel("x_1 = \Delta x_1 + x_{o1} (mol/l)")
    subplot(2,2,4), ylabel("x_2 = \Delta x_2 + x_{o2} (^o C)")
    
    %[~, xout] = simulate(exo_cstr.model, exp_param.t, exp_param.uout, exp_param.x_0);
    
    figure(2)
    subplot(1,2,1)
    plot(exp_param.t, exp_param.r(1,:), 'linestyle', '--', 'color', [0.3 0.3 0.3]); hold on;
    plot(exp_param.t, exp_param.xout(2,:), 'linewidth', 1.5, 'color', cpal(9+i,:)); hold on;     
%     s = scatter(exp_param.t, exp_param.xout(2,:), 'x', 'MarkerEdgeColor', cpal(8,:)); hold on;     
%     s.MarkerFaceAlpha = 0.4;
%     s.MarkerEdgeAlpha = 0.4;
    
    subplot(1,2,2)
    plot(exp_param.t, exp_param.r(2,:), 'linestyle', '--', 'color', [0.3 0.3 0.3]); hold on;
    plot(exp_param.t, exp_param.xout(3,:), 'linewidth', 1.5, 'color', cpal(9+i,:)); hold on;     
%     s = scatter(exp_param.t, exp_param.xout(3,:), 'x', 'MarkerEdgeColor', cpal(8,:)); hold on;     
%     s.MarkerFaceAlpha = 0.4;
%     s.MarkerEdgeAlpha = 0.4;
    end

catch err
    fprintf("Erro!\n");
    rethrow(err)
    
end

% - Exporting the Visualization to an Image
figure(1)
subplot(2,2,1), hold off
subplot(2,2,2), hold off
subplot(2,2,3), hold off
subplot(2,2,4), hold off

timeNow = datetime('now', 'TimeZone', 'local', 'Format', 'dMMMy_HHmmssZ');
figname = "report_experiments/figs/exoSim_control_" + char(timeNow);
fig = gcf; fig.PaperPositionMode = 'auto'; 
print('-bestfit', figname, '-dpdf', '-r300')
system("pdfcrop " + figname + ".pdf " + figname + ".pdf");

figure(2)
subplot(1,2,1) 
ylabel("Concentration (mol/l)")
xlabel("Time (1/hr)")
hold off

subplot(1,2,2)
ylabel("Temperature (^o C)")
xlabel("Time (1/hr)")
hold off

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
foldername = "data/simulations/lqg/";
folder = dir(foldername);
files = {folder.name}; files = files(3:end);
for name = files
    clear exp_param exo_sim simulation
    load(foldername+name)
    
    if(exist("exp_param", "var"))
        for i = 68:-1:1
            if(size(exp_param{i}.xout,1) > 0)
                simulation = exp_param{i};
                save(foldername+name, 'simulation')
                break
            end
        end
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