%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simplified Continuous Stirred Reactor Tank
% 
% Author: Otacilio Bezerra Leite Neto
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd /home/minhotmog/Dropbox/Research/TCC/Codes/exothermal_cstr/
clc; clear all; close all;

% Some colors
cols = [209 17 65;    % 1 - Metro Red
        0 177 89;     % 2 - Metro Green
        0 174 219;    % 3 - Metro Blue
        243 119 53;   % 4 - Metro Orange
        255 196 37;   % 5 - Metro Yellow
        ]/255;  

    
    
%% MODEL PARAMETERS %%
alpha = 30.8285 / 60; beta = 86.688 / 60; gamma = 0.1; delta = 3.556e-4;
K10 = 1.287e12 / 60; K20 = 9.043e6 / 60;
E1 = 9758.3; E2 = 8560;
dH_AB = 4.2; dH_BC = -11; dH_AD = -41.85;
c_in = 5100; T_in = 104.9;

%% REACTION AUXILIARY EQUATIONS %%
K1 = @(T)         K10 * exp( -E1 / (T + 273.15) );
K2 = @(T)         K20 * exp( -E2 / (T + 273.15) );
h  = @(cA, cB, T) -delta * (K1(T)*(cA * dH_AB + cB * dH_BC) + K2(T) * cA^2 * dH_AD);

%% SYSTEM DYNAMICS %%
d_cA = @(t,U,X) - K1(X(3)) * X(1) - K2(X(3)) * X(1)^2 + (c_in - X(1)) * U(1);
d_cB = @(t,U,X)   K1(X(3)) * (X(1) - X(2)) - X(2) * U(1);
d_cC = @(t,U,X)   K1(X(3)) * X(2) - X(5) * U(1);
d_cD = @(t,U,X)   1/2 * K2(X(3)) * X(1)^2 - X(6) * U(1);

d_T  = @(t,U,X)   h(X(1), X(2), X(3)) + alpha * (X(4) - X(3)) + (T_in - X(3)) * U(1);
d_Tc = @(t,U,X)   beta * (X(3) - X(4)) + gamma * U(2);

d_X  = @(t,U,X) [d_cA(t,U,X) d_cB(t,U,X) d_T(t,U,X) d_Tc(t,U,X) d_cC(t,U,X) d_cD(t,U,X)];

%% Simulation (Using Euler Method) %%


%% STATIONARY POINTS %%
% Inputs Stationary Spaces
U_ss = [linspace(5/60, 35/60, 100)' linspace(-8.5/60, 0, 100)'];

% Stationary Space Numerical Calculation
% X_ss = zeros(100,100,4);
% prev = [0 0 0 0];
% for i = 1:1:100
%     for j = 1:1:100
%         F = @(X) [-K1(X(3))*X(1) - K2(X(3))*X(1)^2 + (c_in - X(1)) * U_ss(i,1);
%                    K1(X(3))*(X(1) - X(2)) - X(2)*U_ss(i,1);
%                    h(X(1),X(2),X(3)) + alpha*(X(4) - X(3)) + (T_in - X(3))* U_ss(i,1);
%                    beta*(X(3) - X(4)) + gamma*U_ss(j,2)];      
%         
%         X_ss(i,j,:) = fsolve(F, prev);
%         prev = Y(i,j,:);
%     end
% end

lStruct = load('stationaryPoints.mat');
X_ss = lStruct.Y; 

% Visualization of the Stationary Space
figure;
subplot(1,2,1), mesh(U_ss(:,1), U_ss(:,2), X_ss(:,:,2)')
title("Stationary Points - C_B")
xlabel("Input Flow-rate (m^3/min)"), ylabel("Cooling Capacity (MJ/h)"), zlabel("Outflow Concentration (mol/m^3)")

subplot(1,2,2), mesh(U_ss(:,1), U_ss(:,2), X_ss(:,:,3)')
title("Stationary Points - T")
xlabel("Input Flow-rate (m^3/min)"), ylabel("Cooling Capacity (MJ/h)"), zlabel("Tank Temperature (K)")

% Operating Space
lStruct = load('opPoints.mat');
U_ss = lStruct.U_ss;
X_ss = lStruct.X_ss;

figure;
subplot(1,2,1), mesh(U_ss(:,1), U_ss(:,2), X_ss(:,:,2)', H)
title("Operation Points - C_B")
xlabel("Input Flow-rate (m^3/min)"), ylabel("Cooling Capacity (MJ/h)"), zlabel("Outflow Concentration (mol/m^3)")

subplot(1,2,2), mesh(U_ss(:,1), U_ss(:,2), X_ss(:,:,3)', H)
title("Operation Points - T")
xlabel("Input Flow-rate (m^3/min)"), ylabel("Cooling Capacity (MJ/h)"), zlabel("Tank Temperature (K)")

%% LINEAR STATE-SPACE REPRESENTATION %%
% Auxiliary Derivatives
dK_1  = @(Te)         -(K10 * E1 * exp(-E1 / (Te + 273.12))) / (Te +273.12)^2; 
dK_2  = @(Te)         -(K20 * E2 * exp(-E2 / (Te + 273.12))) / (Te +273.12)^2;
dh_dt = @(cAe,cBe,Te) -gamma*(dK_1(Te)*(cAe*dH_AB + cBe*dH_BC) + dK_2(Te)*cAe^2*dH_AD);

% Linearized Matrices
A_e = @(U_e,X_e) [-U_e(1)-K1(X_e(3))-2*K2(X_e(3))*X_e(1)         0                          dK_1(X_e(3))*X_e(1)-dK_2(X_e(3))*X_e(1)^2    0     ;
                   K1(X_e(3))                                   -U_e(1)-K1(X_e(3))          dK_1(X_e(3))*(X_e(1) - X_e(2))               0     ;  
                  -gamma*(K1(X_e(3))*dH_AB+2*K2(X_e(3))*dH_AB)  -gamma*(K1(X_e(3))*dH_BC)  -U_e(1)-alpha+dh_dt(X_e(1),X_e(2),X_e(3))     alpha ;
                   0                                             0                          beta                                        -beta  ];

B_e = @(U_e,X_e) [ c_in - X_e(1)      0;
                    -X_e(2)          0;
                  T_in - X_e(3)      0;
                       0         gamma];

C = [0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0];
D = [0 0; 0 0; 0 0; 0 0];

%% SYSTEM RESPONSE MODES %%
% Calculating the modes
sys_modes = cell(size(X_ss));
sys_modes{1,1} = 1;

lambda = zeros(size(X_ss));

for i = 1:size(X_ss,1)
    for j = 1:size(X_ss,2)
        lambda(i,j,:) = eig(A_e([U_ss(i,1) U_ss(j,2)], X_ss(i,j,:)));

        sys_modes{i,j,1} = eval(sprintf('@(t) exp(%d * t)', lambda(i,j,1)));
        sys_modes{i,j,2} = eval(sprintf('@(t) exp(%d * t) .* cos(%d * t)', real(lambda(i,j,2)), imag(lambda(i,j,2))));
        sys_modes{i,j,3} = sys_modes{i,j,2};
        sys_modes{i,j,4} = eval(sprintf('@(t) exp(%d * t)', lambda(i,j,4)));
    end
end

% Visualizing the modes
t = 0:0.1:8;
figure;
for i = 1:size(X_ss,1) 
    for j = 1:size(X_ss,2)
        tic
        subplot(1,3,1), plot(t, sys_modes{i,j,1}(t), 'color', H(i,j,:)), title("1º Mode"), hold on
        subplot(1,3,2), plot(t, sys_modes{i,j,2}(t), 'color', H(i,j,:)), title("2º and 3º Mode"), hold on
        subplot(1,3,3), plot(t, sys_modes{i,j,4}(t), 'color', H(i,j,:)), title("4º Mode"), hold on
        drawnow
        toc
    end
end

% Matrix Exponential Visualization
syms t
ts = 0:0.1:5;
e_At = zeros(size(X_ss));

figure;
for i = 1:size(X_ss,1) 
    for j = 1:size(X_ss,2)
        [V, J] = jordan(A_e([U_ss(i,1) U_ss(j,2)], X_ss(i,j,:)));
        e_At = V*expm(J*t)*pinv(V);

        tic
        for k = 1:4
            subplot(4,4,k), plot(ts, double(subs(e_At(1,k),ts)), 'color', H(i,j,:)), hold on
            subplot(4,4,4+k), plot(ts, double(subs(e_At(2,k),ts)), 'color', H(i,j,:)), hold on
            subplot(4,4,8+k), plot(ts, double(subs(e_At(3,k),ts)), 'color', H(i,j,:)), hold on
            subplot(4,4,12+k), plot(ts, double(subs(e_At(4,k),ts)), 'color', H(i,j,:)), hold on
        end
        drawnow
        toc
    end
end

save('modeVisualization.mat', 'sys_modes', 'lambda', 'e_At')

%% LINEAR STATE-SPACE SIMULATION %%
% Selects Linearized Model
idx = 1;
A = A_e(U_ss(idx,:), X_ss(idx,idx,:)); B = B_e(U_ss(idx,:), X_ss(idx,idx,:));
cA_ss = X_ss(idx,idx,1); cB_ss = X_ss(idx,idx,2); T_ss = X_ss(idx,idx,3); Tc_ss = X_ss(idx,idx,4);
U1_ss = U_ss(idx,1); U2_ss = U_ss(idx,2);

d_X = @(t,U,X) [d_cA(t,U,X) d_cB(t,U,X) d_T(t,U,X) d_Tc(t,U,X)];

% Simulation Parameters
t  = 0:0.1:30;                                                 % Time 
y0 = [cA_ss cB_ss T_ss Tc_ss];                                 % Initial Conditions
U  = [(sin(t)*0.05 + U1_ss)' (sin(t)*0.05 + U2_ss)'];          % Input Signal

% Simulation of the Outputs (Linear)
[~, ~, yout1] = sim('exothermalCSTR_LinearSS', t, [], [t', U]);

% Simulation of the Outputs (Non-Linear)
yout2 = zeros(size(t,1), size(y0,2));
yout2(1,:) = y0;
for i = 2:size(t,2)
    [~, y_aux] = odeSolver(d_X, t(i-1:i), U(i-1,:), yout2(i-1,:), 100); 
    yout2(i,:) = y_aux(end,:);
end

% Visualization
figure;
subplot(2,3,1), plot(t,U(:,1)), title("Input Signal")
xlabel("Time (min)"), ylabel("Input Flow-rate (m^3/min)")
grid()

subplot(2,3,2), plot(t,yout2(:, 1:2)), title("Output Signals (Linear)")
legend("C_A", "C_B"),
xlabel("Time (min)"), ylabel("Outflow Concentration (mol/m^3)")
grid()

subplot(2,3,3), plot(t,yout1(:, 1:2)), title("Output Signals (Linear)")
legend("C_A", "C_B"),
xlabel("Time (min)"), ylabel("Outflow Concentration (mol/m^3)")
grid()

subplot(2,3,4), plot(t,U(:,2)), title("Input Signal")
xlabel("Time (min)"), ylabel("Cooling Capacity (MJ/min)")
grid()

subplot(2,3,5), plot(t,yout2(:, 3:4)), title("Output Signals")
legend("T", "T_C"), xlabel("Time (min)"), ylabel("Temperatures (K)")
grid()

subplot(2,3,6), plot(t,yout1(:, 3:4)), title("Output Signals (Linear)")
legend("T", "T_C"), xlabel("Time (min)"), ylabel("Temperatures (K)")
grid()

%% LINEAR STABILITY ANALYSIS %%
figure;
for i = 1:size(X_ss,1)
    for j = 1:size(X_ss,2)
        tic
        A = A_e([U_ss(j,1) U_ss(j,2)], X_ss(i,j,:));
        plot(eig(A), 'x', 'color', H(i,j,:), 'markersize', 10); hold on
        drawnow
        toc
    end
end
title("Pole-Zero Map")
xlabel("Real Axis (secons^{-1})")
ylabel("Imaginary Axis (secons^{-1})")
sgrid()

figure;
for i = 1:size(X_ss,1)
    for j = 1:size(X_ss,2)
        tic
        A = A_e([U_ss(j,1) U_ss(j,2)], X_ss(i,j,:));
        B = B_e([U_ss(j,1) U_ss(j,2)], X_ss(i,j,:));
        
        [NUM, DEN] = ss2tf(A,B,C,D,1);
        G = [tf(NUM(2,:), DEN); tf(NUM(3,:), DEN)];

        bodeplot(G,'b');
        
        lineHandle = findobj(gcf,'Type','line','-and','Color','b');
        set(lineHandle,'Color',H(i,j,:));

        hold on
        drawnow
        toc
        i*j
    end
end
grid()
saveas(gcf, '../Presentations/September_2018/imgs/exoCSTR_simu09', 'epsc')
close all

figure;
for i = 1:size(X_ss,1)
    for j = 1:size(X_ss,2)
        tic
        A = A_e([U_ss(j,1) U_ss(j,2)], X_ss(i,j,:));
        B = B_e([U_ss(j,1) U_ss(j,2)], X_ss(i,j,:));
        
        [NUM, DEN] = ss2tf(A,B,C,D,1);
        G = [tf(NUM(2,:), DEN); tf(NUM(3,:), DEN)];

        nyquistplot(G,'b');
        
        lineHandle = findobj(gcf,'Type','line','-and','Color','b');
        set(lineHandle,'Color',H(i,j,:));

        hold on
        drawnow
        toc
        i*j
    end
end
grid()
saveas(gcf, '../Presentations/September_2018/imgs/exoCSTR_simu10', 'epsc')
close all

%% CONTROLABILITY AND OBSERVABILITY %%
idx = 1;
A = A_e(U_ss(idx,:), X_ss(idx,idx,:)); B = B_e(U_ss(idx,:), X_ss(idx,idx,:));
C = [0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0];

CO = ctrb(A,B);
OB = obsv(A,C);