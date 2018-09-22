%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simplified Continuous Stirred Reactor Tank
% 
% Author: Otacilio Bezerra Leite Neto
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd /home/minhotmog/Dropbox/Research/TCC/Codes/isothermal_cstr/
clc; clear all; close all;

% Some colors
cols = [209 17 65;    % 1 - Metro Red
        0 177 89;     % 2 - Metro Green
        0 174 219;    % 3 - Metro Blue
        243 119 53;   % 4 - Metro Orange
        255 196 37;   % 5 - Metro Yellow
        ]/255;  
 
%% MODEL PARAMETERS %%
K1 = 5/6; K2 = 20/3; K3 = 1/6;
cAf = 10; V = 1;

%% CONCENTRATION DYNAMICS %%
d_cA = @(t,U,X)  U(1)/V * (cAf - X(1)) + (-K1 * X(1) - K3 * X(1)^2);
d_cB = @(t,U,X) -U(1)/V * X(2)         + (K1 * X(1) - K2 * X(2));
d_cC = @(t,U,X) -U(1)/V * X(3)         + (K2 * X(2));
d_cD = @(t,U,X) -U(1)/V * X(4)         + (1/2 * K3 * X(1)^2);

d_C = @(t,U,X) [d_cA(t,U,X) d_cB(t,U,X) d_cC(t,U,X) d_cD(t,U,X)];

%% Simulation (Using Euler Method) %%
% Simulation Parameters
t = (0:0.1:10)';                                    % Time
y0 = [0 0 0 0];                                     % Initial Conditions
U = [ones(15,1)*0; ones(40,1); ones(numel(t)-55,1)*0] ;   % Input Signal

% Simulation of the Outputs
y = zeros(size(t,1), size(y0,2));
y(1,:) = y0;
for i = 2:size(t,1)
    [~, y_aux] = odeSolver(d_C, t(i-1:i), U(i-1), y(i-1,:), 100); 
    y(i,:) = y_aux(end,:);
end
    
% Visualization
figure;
subplot(1,2,1), plot(t,U, 'linewidth', 1.5), title("Input Signal")
xlabel("Time (min)"), ylabel("Input Flow-rate (m^3/min)")
grid()

subplot(1,2,2), p = plot(t,y); title("Output Signals")
p(1).LineWidth = 1.5; p(2).LineWidth = 1.5;
p(3).LineStyle='--'; p(4).LineStyle='--';
legend("C_A", "C_B", "C_C", "C_D") 
xlabel("Time (min)"), ylabel("Outflow Concentration (mol/m^3)")
grid()

%% Simulation (Using Simulink) %%
% Simulation Parameters
t = 0:0.1:10;                                       % Time 
y0 = [0 0 0 0];                                     % Initial Conditions
U = ones(numel(t),1) * 3.03;                     % Input Signal

% Simulation of the Outputs
[~, ~, yout] = sim('isothermalCSTR', t, [], [t', U]);
    
% Visualization
figure;
subplot(1,2,1), title("Input Signal"), plot(t,U)
xlabel("Time (min)"), ylabel("Input Flow-rate (m^3/min)")
grid()

subplot(1,2,2), title("Output Signals"), p = plot(t,y);
p(1).LineWidth = 1.5; p(2).LineWidth = 1.5;
p(3).LineStyle='--'; p(4).LineStyle='--';
legend("C_A", "C_B", "C_C", "C_D") 
xlabel("Time (min)"), ylabel("Outflow Concentration (mol/m^3)")
grid()

%% OPERATING POINTS %%
% Stationary Equations
cA_ss = @(U) (-(K1 + U) + sqrt((K1 + U).^2 + 4 * K3 * U * cAf)) / (2 * K3);
cB_ss = @(U) (K1 .* cA_ss(U)) ./ (K2 + U);
cC_ss = @(U) (K2 .* cB_ss(U)) ./ U;
cD_ss = @(U) (1/2 * K3 * cA_ss(U).^2) ./ U;

% Input Stationary Space
U_ss = linspace(0, 20, 100);

% Output Stationary Space
Y_ss = cB_ss(U_ss);

% Visualization of Operating Points
figure; 
plot(U_ss, Y_ss,'linewidth',1.5), grid()
title("Stationary Points - C_B")
xlabel("Input Flow-rate (m^3/min)"), ylabel("Outflow Concentration (mol/m^3)")

% Getting the Operating Space
[~, I] = max(Y_ss);
U_ss = U_ss(I-10:I+10)';
X_ss = [cA_ss(U_ss) cB_ss(U_ss)];

blue_S = [14,104,206]/255;
blue_E = [244,71,71]/255;
colors_p = [linspace(blue_S(1),blue_E(1),21)', linspace(blue_S(2),blue_E(2),21)', linspace(blue_S(3),blue_E(3),21)'];

figure; 
hold on
for n = 1:20
    plot(U_ss([n n+1]), X_ss([n n+1],2), 'color', colors_p(n,:), 'linewidth', 1.5);
end
title("Stationary Points - C_B")
xlabel("Input Flow-rate (m^3/min)"), ylabel("Outflow Concentration (mol/m^3)")
grid()

%% LINEAR STATE-SPACE REPRESENTATION %%
A_e = @(U_e,X_e) [-U_e(1) - K1 - 2*K3*X_e(:,1)          0     ;
                         K1                    -U_e(1) - K2];

B_e   = @(U_e,X_e) [ cAf - X_e(:,1);
                     -X_e(:,2)   ];

C = eye(2);
D = [0; 0;];

%% SYSTEM RESPONSE MODES %%
% Calculating the modes
sys_modes = cell(size(X_ss));
sys_modes{1,1} = 1;
for i = 1:size(X_ss,1)
    lambda = eig(A_e(U_ss(i), X_ss(i,:)));
    
    sys_modes{i,1} = eval(sprintf('@(t) exp(%d * t)', lambda(1)));
    sys_modes{i,2} = eval(sprintf('@(t) exp(%d * t)', lambda(2)));
end

% Visualizing the modes
t = 0:0.1:5;

figure;
for i = 1:size(X_ss,1)
    subplot(1,2,1), plot(t, sys_modes{i,1}(t), 'color', colors_p(i,:)), title("First Mode"), hold on
    subplot(1,2,2), plot(t, sys_modes{i,2}(t), 'color', colors_p(i,:)), title("Second Mode"), hold on
end

% Matrix Exponential Visualization
syms t
ts = 0:0.1:5;

figure;
for i = 1:size(X_ss,1)
    e_At = expm(A_e(U_ss(i), X_ss(i,:))*t);
    
    subplot(2,2,1), plot(ts, double(subs(e_At(1,1),ts)), 'color', colors_p(i,:)), hold on
    subplot(2,2,2), plot(ts, double(subs(e_At(1,2),ts)), 'color', colors_p(i,:)), hold on
    subplot(2,2,3), plot(ts, double(subs(e_At(2,1),ts)), 'color', colors_p(i,:)), hold on
    subplot(2,2,4), plot(ts, double(subs(e_At(2,2),ts)), 'color', colors_p(i,:)), hold on
end

%% LINEAR STATE-SPACE SIMULATION %%
% Selects Linearized Model
idx = 1;
A = A_e(U_ss(idx,:), X_ss(idx,:)); B = B_e(U_ss(idx,:), X_ss(idx,:));
U_ss = U_ss(idx,1);
cA_e = X_ss(idx,1); cB_e = X_ss(idx,2);

% Simulation Parameters
t  = 0:0.1:30;                             % Time 
y0 = [cA_e cB_e cC_e cD_e];            % Initial Conditions
U  = (sin(t) + U_ss)';                     % Input Signal

% Simulation of the Outputs (Linear)
[~, ~, yout1] = sim('isothermalCSTR_LinearSS', t, [], [t', U]);

% Simulation of the Outputs (Non-Linear)
[~, ~, yout2] = sim('isothermalCSTR', t, [], [t', U]);

% Visualization
figure;
subplot(1,3,1), title("Input Signal"), plot(t,U)
xlabel("Time (min)"), ylabel("Input Flow-rate (m^3/min)")
grid()

subplot(1,3,2), title("Output Signals (Linearized)"), plot(t,yout1);
legend("C_A", "C_B") 
xlabel("Time (min)"), ylabel("Outflow Concentration (mol/m^3)")
grid()

subplot(1,3,3), title("Output Signals (Non-Linear)"), p = plot(t,yout2);
p(1).LineWidth = 1.5; p(2).LineWidth = 1.5;
p(3).LineStyle='--'; p(4).LineStyle='--';
legend("C_A", "C_B", "C_C", "C_D") 
xlabel("Time (min)"), ylabel("Outflow Concentration (mol/m^3)")
grid()

%% LINEAR STABILITY ANALYSIS %%
blue_S = [14,104,206]/255;
blue_E = [244,71,71]/255;
colors_p = [linspace(blue_S(1),blue_E(1),22)', linspace(blue_S(2),blue_E(2),22)', linspace(blue_S(3),blue_E(3),22)'];

figure;
for j = 1:size(X_ss,1)
    A = A_e(U_ss(j), X_ss(j,:));
    B = B_e(U_ss(j), X_ss(j,:));
    
    [NUM, DEN] = ss2tf(A,B,C,D);
    G = [tf(NUM(1,:), DEN); tf(NUM(2,:), DEN)];
    
    pzmap(G); hold on   
end

a = findobj(gca,'type','line');
for i = 1:length(a)
    set(a(i), 'markersize',10)
    set(a(i), 'linewidth',2)
    set(a(i), 'color', colors_p(floor(i/2)+1,:))
end
grid()

figure;
for j = 1:size(X_ss,1)
    A = A_e(U_ss(j), X_ss(j,:));
    B = B_e(U_ss(j), X_ss(j,:));
    
    [NUM, DEN] = ss2tf(A,B,C,D);
    G = [tf(NUM(1,:), DEN); tf(NUM(2,:), DEN)];
    
    nyquistplot(G,'b');

    lineHandle = findobj(gcf,'Type','line','-and','Color','b');
    set(lineHandle,'Color',colors_p(j,:));

    hold on  
end
grid()

%% CONTROLABILITY AND OBSERVABILITY %%
idx = 1;
A = A_e(U_ss(idx,:), X_ss(idx,:)); B = B_e(U_ss(idx,:), X_ss(idx,:));
C = [0 0; 0 1];

CO = ctrb(A,B);
OB = obsv(A,C);