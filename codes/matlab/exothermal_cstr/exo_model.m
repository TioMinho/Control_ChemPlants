%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exothermal Continuous Stirred Reactor Tank
% 
% Author: Otacilio Bezerra Leite Neto
%
%   This script defines and saves a model for an Exothermal CSTR.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%s\n', "Loading the Exothermal CSTR model...")
vars = who;
%% MODEL PARAMETERS %%
alpha = 10 / 60; beta = 100 / 60; gamma = 0.3; delta = 3.556e-4;
K10 = 1.287e12 / 60; K20 = 9.043e6 / 60;
E1 = 9758.3; E2 = 9000;
dH_AB = 4.2; dH_BC = -30; dH_AD = -41.85;
c_in = 5100; T_in = 105;

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

d_X  = @(t,U,X) [d_cA(t,U,X) d_cB(t,U,X) d_T(t,U,X) d_Tc(t,U,X)];
sysVar = @(t,U,X) [d_cA(t,U,X) d_cB(t,U,X) d_T(t,U,X) d_Tc(t,U,X) d_cC(t,U,X) d_cD(t,U,X)];

%% OPERATING POINTS %%
% Stationary Space Numerical Calculation
U_ss = [linspace(5/60, 35/60, 25)' linspace(-8.5/60, 0, 25)'];
X_ss = zeros(25,25,4);
prev = [0 0 0 0];
for i = 1:1:25
    for j = 1:1:25
        F = @(X) [-K1(X(3))*X(1) - K2(X(3))*X(1)^2 + (c_in - X(1)) * U_ss(i,1);
                   K1(X(3))*(X(1) - X(2)) - X(2)*U_ss(i,1);
                   h(X(1),X(2),X(3)) + alpha*(X(4) - X(3)) + (T_in - X(3))* U_ss(i,1);
                   beta*(X(3) - X(4)) + gamma*U_ss(j,2)];      
        
        X_ss(i,j,:) = fsolve(F, prev);
        prev = X_ss(i,j,:);
    end
end

% Operating Space
% lStruct = load('../data/exoCSTR_opPoints.mat');
% U_ss = lStruct.U_ss;
% X_ss = lStruct.X_ss;
% numOp = 25*25;

figure(3);
subplot(1,2,1), surf(U_ss(:,1), U_ss(:,2), X_ss(:,:,2)', ...
                                    ccmap.area(:,:,1:3), 'AlphaData', ccmap.area(:,:,4), 'FaceAlpha', 'flat')
title("Operation Points - C_B")
xlabel("Input Flow-rate (m^3/min)"), ylabel("Cooling Capacity (MJ/h)"), zlabel("Outflow Concentration (mol/m^3)")

subplot(1,2,2), surf(U_ss(:,1), U_ss(:,2), X_ss(:,:,3)', ...
                                    ccmap.area(:,:,1:3), 'AlphaData', ccmap.area(:,:,4), 'FaceAlpha', 'flat')
title("Operation Points - T")
xlabel("Input Flow-rate (m^3/min)"), ylabel("Cooling Capacity (MJ/h)"), zlabel("Tank Temperature (K)")


%% LINEAR STATE-SPACE REPRESENTATION %%
% Auxiliary Derivatives
dK_1  = @(Te)         -(K10 * E1 * exp(-E1 / (Te + 273.12))) / (Te +273.12).^2; 
dK_2  = @(Te)         -(K20 * E2 * exp(-E2 / (Te + 273.12))) / (Te +273.12).^2;
dh_dt = @(cAe,cBe,Te) -gamma*(dK_1(Te)*(cAe*dH_AB + cBe*dH_BC) + dK_2(Te)*cAe^2*dH_AD);

posX = @(i) mod(i,25)+(mod(i,25) == 0)*25;
posY = @(i) ceil(i/25);

% Linearized Matrices
A_e = @(op) [-U_ss(posX(op), 1)-K1(X_ss(posX(op), posY(op),3))-2*K2(X_ss(posX(op), posY(op),3))*X_ss(posX(op), posY(op),1)         0                          dK_1(X_ss(posX(op), posY(op),3))*X_ss(posX(op), posY(op),1)-dK_2(X_ss(posX(op), posY(op),3))*X_ss(posX(op), posY(op),1)^2    0     ;
                   K1(X_ss(posX(op), posY(op),3))                                   -U_ss(posX(op),1)-K1(X_ss(posX(op), posY(op),3))          dK_1(X_ss(posX(op), posY(op),3))*(X_ss(posX(op), posY(op),1) - X_ss(posX(op), posY(op),2))               0     ;  
                  -gamma*(K1(X_ss(posX(op), posY(op),3))*dH_AB+2*K2(X_ss(posX(op), posY(op),3))*dH_AB)  -gamma*(K1(X_ss(posX(op), posY(op),3))*dH_BC)  -U_ss(posX(op),1)-alpha+dh_dt(X_ss(posX(op), posY(op),1),X_ss(posX(op), posY(op),2),X_ss(posX(op), posY(op),3))     alpha ;
                   0                                             0                          beta                                        -beta  ];

B_e = @(op) [ c_in - X_ss(posX(op), posY(op),1)      0;
                    -X_ss(posX(op), posY(op),2)          0;
                  T_in - X_ss(posX(op), posY(op),3)      0;
                       0         gamma];

C = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
D = [0 0; 0 0; 0 0; 0 0];

%% SYSTEM POLES FROM THE CHARACTERISTIC POLYNOMIAL %%
lambda = zeros(numOp, 4);
for i = 1:numOp
    lambda(i,:) = eig(A_e(i));
end

%% SYSTEM MODES %%
modes = @(t) exp(real(lambda)*t) .* cos(imag(lambda)*t);
modes = sym(modes);

%% SYSTEM TRANSFER MATRICES %%
syms t
e_At = sym(t * ones(4, 4, numOp));
for i = 1:numOp
    [V, J] = jordan(A_e(i));
    e_At(:,:,i) = V*expm(J*t)*pinv(V);
end

%% SAVES THE MODEL %%
newVars = setdiff(who, vars);

exo_cstr.param         = struct('alpha', alpha, 'beta', beta, 'gamma', gamma, 'delta', delta, 'K10', K10, 'K20', K20, ...
                                        'E1', E1, 'E2', E2, 'dH_AB', dH_AB, 'dH_BC', dH_BC, 'dH_AD', dH_AD, 'c_in', c_in, 'T_in', T_in);
exo_cstr.model         = d_X;
exo_cstr.sysVar         = sysVar;
exo_cstr.oper            = struct('U', U_ss, 'X', X_ss, 'size', size(X_ss, 1)*size(X_ss, 2));
exo_cstr.ss_model    = struct('A', A_e, 'B', B_e, 'C', C, 'D', D);
exo_cstr.poles           = lambda;
exo_cstr.modes        = modes;
exo_cstr.trf_matrix    = e_At;
exo_cstr.sizeX           = size(C, 1);
exo_cstr.sizeU          = size(D, 2);
exo_cstr.sizeY           = sum(sum(C, 2) ~= 0);

save('../data/exo_cstr_model.mat', 'exo_cstr')

% Clean up the mess
clear(newVars{:})
clear newVars
