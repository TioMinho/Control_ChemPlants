%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simplified Continuous Stirred Reactor Tank
% 
% Author: Otacilio Bezerra Leite Neto
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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