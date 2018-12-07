%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Isothermal Continuous Stirred Reactor Tank
% 
% Author: Otacilio Bezerra Leite Neto
%
%   This script defines and saves a model for an Isothermal CSTR.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%s\n', "Loading the Isothermal CSTR model...")
vars = who;
%% MODEL PARAMETERS %%
K1 = 5/6; K2 = 5/3; K3 = 1/6;
cAf = 10; V = 1;

%% SYSTEM DYNAMICS %%
d_cA = @(t,U,X)  U(1) * (cAf - X(:,1)) + (-K1 * X(:,1) - K3 * X(:,1)^2);
d_cB = @(t,U,X) -U(1) * X(:,2)            + (K1 * X(:,1) - K2 * X(:,2));
d_cC = @(t,U,X) -U(1) * X(:,3)            + (K2 * X(:,2));
d_cD = @(t,U,X) -U(1) * X(:,4)            + (1/2 * K3 * X(:,1)^2);

d_X   = @(t,U,X) [d_cA(t,U,X) d_cB(t,U,X)];
sysVar = @(t,U,X) [d_cA(t,U,X) d_cB(t,U,X) d_cC(t,U,X) d_cD(t,U,X)];

%% OPERATING POINTS %%
% Stationary Equations
cA_ss = @(U) (-(K1 + U) + sqrt((K1 + U).^2 + 4 * K3 * U * cAf)) / (2 * K3);
cB_ss = @(U) (K1 .* cA_ss(U)) ./ (K2 + U);
cC_ss = @(U) (K2 .* cB_ss(U)) ./ U;
cD_ss = @(U) (1/2 * K3 * cA_ss(U).^2) ./ U;

% Input Stationary Space
U_ss = linspace(0, 20, 200);

% Output Stationary Space
Y_ss = cB_ss(U_ss);

% % Visualization of the Stationary Space
% figure; 
% plot(U_ss, Y_ss,'linewidth',1.5), grid()
% title("Stationary Points - C_B")
% xlabel("Input Flow-rate (m^3/min)"), ylabel("Outflow Concentration (mol/m^3)")

% Getting the Operating Space
I = find(Y_ss > 2/3 * max(Y_ss));
U_ss = U_ss(I)';
X_ss = [cA_ss(U_ss) cB_ss(U_ss)];

% % Visualization of the Operation Space
% figure; 
% plot(U_ss, X_ss(:,2),'linewidth',1.5), grid()
% title("Stationary Points - C_B")
% xlabel("Input Flow-rate (m^3/min)"), ylabel("Outflow Concentration (mol/m^3)")

%% LINEAR STATE-SPACE REPRESENTATION %%
A_e = @(op)      [-U_ss(op,1) - K1 - 2*K3*X_ss(op,1)                  0      ;
                                                 K1                             -U_ss(op,1) - K2];

B_e   = @(op)    [ cAf - X_ss(op,1);
                                      -X_ss(op,2)   ];

C = [1 0; 0 1];
D = [0; 0];

%% SYSTEM POLES FROM THE CHARACTERISTIC POLYNOMIAL %%
lambda = zeros(size(X_ss));
for i = 1:size(X_ss,1)
    lambda(i,:) = eig(A_e(i));
end

%% SYSTEM MODES %%
modes = @(t) exp(lambda*t);
modes = sym(modes);

%% SYSTEM TRANSFER MATRICES %%
syms t
e_At = sym(t * ones(2, 2, size(X_ss,1)));
for i = 1:size(X_ss,1)
    e_At(:,:,i) = expm(A_e(i) * t);
end

%% CONSTRUCTS THE MODEL %%
newVars = setdiff(who, vars);

iso_cstr.param        = struct('K1', K1, 'K2', K2, 'K3', K3, 'cA_f', cAf);
iso_cstr.model        = d_X;
iso_cstr.sysVar       = sysVar;
iso_cstr.oper         = struct('U', U_ss, 'X', X_ss, 'size', size(X_ss,1));
iso_cstr.ss_model     = struct('A', A_e, 'B', B_e, 'C', C, 'D', D);
iso_cstr.poles        = lambda;
iso_cstr.modes        = modes;
iso_cstr.trf_matrix   = e_At;
iso_cstr.sizeX        = size(C, 1);
iso_cstr.sizeU        = size(D, 2);
iso_cstr.sizeY        = sum(sum(C, 2) ~= 0);

% Clean up the mess
clear(newVars{:})
clear newVars