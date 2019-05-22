%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamical Modelling and Intelligent Control of Chemical Reactive Systems
% 
% Author: Otacilio Bezerra Leite Neto
%
%   experiment_parameters.m
%       This script contains the a series of parameters used for simulation
%       the system several times.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_experiments = 69;

exp_param = cell(n_experiments);

%1 
i = 1;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/4).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)].*[0.9;0.99] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([10, 10, 10, 10, 1e14, 1e14]);
R = diag([10 10]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe;
exp_param{i}.type = 'lqri'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%2 
i = 2;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/4).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)].*[0.9;0.99] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([1, 1, 1, 1, 1e14, 1e14]);
R = diag([1 1]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe;
exp_param{i}.type = 'lqri'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%3 
i = 3;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/4).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)].*[0.9;0.99] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([1, 1, 1, 1, 1e15, 1e15]);
R = diag([1 1]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe;
exp_param{i}.type = 'lqri'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%4 
i = 4;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/4).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)].*[0.9;0.99] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([100, 100, 100, 100, 1e14, 1e14]);
R = diag([100 100]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe;
exp_param{i}.type = 'lqri'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%5 
i = 5;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/2).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([1, 1, 1, 1, 1e14, 1e14]);
R = diag([1 1]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe;
exp_param{i}.type = 'lqri'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%6 
i = 6;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/2).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([100, 100, 100, 100, 1e14, 1e14]);
R = diag([100 100]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe;
exp_param{i}.type = 'lqri'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%7 
i = 7;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/2).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([1, 1, 1, 1]);
R = diag([1 1]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe;
exp_param{i}.type = 'lqr'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%8 
i = 8;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/2).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([100, 100, 100, 100]);
R = diag([100 100]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe;
exp_param{i}.type = 'lqr'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%9
i = 9;
t = (0:0.1:23.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([10, 25, 25, 10]);
R = diag([10 10]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1 0.9 0.99 1];
exp_param{i}.type = 'lqr'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%10
i = 10;
t = (0:0.1:23.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([10, 25, 25, 10]);
R = diag([10 10]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [0.8 0.8 0.98 0.98];
exp_param{i}.type = 'lqr'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%11
i = 11;
t = (0:0.1:23.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([1, 1, 1, 1]);
R = diag([1 1]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [0.8 0.8 0.98 0.98];
exp_param{i}.type = 'lqr'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%12
i = 12;
t = (0:0.1:23.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([100, 100, 100, 100]);
R = diag([100 100]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [0.8 0.8 0.98 0.98];
exp_param{i}.type = 'lqr'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%13
i = 13;
t = (0:0.1:23.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([1000, 1000, 100, 100]);
R = diag([100 100]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [0.8 0.8 0.98 0.98];
exp_param{i}.type = 'lqr'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%14
i = 14;
t = (0:0.1:23.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([1000, 1000, 100, 100]);
R = diag([100 100]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1.2 1.2 1.02 1.02];
exp_param{i}.type = 'lqr'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%15
i = 15;
t = (0:0.1:23.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([1, 1, 1, 1]);
R = diag([1 1]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1.2 1.2 1.02 1.02];
exp_param{i}.type = 'lqr'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%16
i = 16;
t = (0:0.1:23.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([10, 10, 10, 10]);
R = diag([10 10]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1.2 1.2 1.02 1.02];
exp_param{i}.type = 'lqr'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%17 
i = 17;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/2).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([1, 1, 1, 1]);
R = diag([1 1]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe;
exp_param{i}.type = 'lqr'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%18 
i = 18;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/2).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([100, 100, 100, 100]);
R = diag([100 100]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe;
exp_param{i}.type = 'lqr'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%19
i = 19;
t = (0:0.1:4.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([10, 25, 25, 10]);
R = diag([10 10]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1 0.9 0.99 1];
exp_param{i}.type = 'lqr'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%20
i = 20;
t = (0:0.1:4.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([10, 25, 25, 10]);
R = diag([10 10]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [0.8 0.8 0.98 0.98];
exp_param{i}.type = 'lqr'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%21
i = 21;
t = (0:0.1:4.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([1, 1, 1, 1]);
R = diag([1 1]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [0.8 0.8 0.98 0.98];
exp_param{i}.type = 'lqr'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%22
i = 22;
t = (0:0.1:4.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([100, 100, 100, 100]);
R = diag([100 100]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [0.8 0.8 0.98 0.98];
exp_param{i}.type = 'lqr'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%23
i = 23;
t = (0:0.1:4.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([1000, 1000, 100, 100]);
R = diag([100 100]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [0.8 0.8 0.98 0.98];
exp_param{i}.type = 'lqr'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%24
i = 24;
t = (0:0.1:4.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([1000, 1000, 100, 100]);
R = diag([100 100]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1.2 1.2 1.02 1.02];
exp_param{i}.type = 'lqr'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%25
i = 25;
t = (0:0.1:4.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([1, 1, 1, 1]);
R = diag([1 1]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1.2 1.2 1.02 1.02];
exp_param{i}.type = 'lqr'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%26
i = 26;
t = (0:0.1:4.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([10, 10, 10, 10]);
R = diag([10 10]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1.2 1.2 1.02 1.02];
exp_param{i}.type = 'lqr'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%27
i = 27;
t = (0:0.1:4.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([10, 10, 10, 10]);
R = diag([10 10]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1.4 1.4 1.04 1.04];
exp_param{i}.type = 'lqr'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%28
i = 28;
t = (0:0.1:4.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([10, 10, 10, 10]);
R = diag([10 10]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [0.6 0.6 0.96 0.96];
exp_param{i}.type = 'lqr'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%29
i = 29;
t = (0:0.1:4.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([10, 10, 10, 10]);
R = diag([10 10]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1 0.6 0.96 1];
exp_param{i}.type = 'lqr'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%30
i = 30;
t = (0:0.1:4.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([10, 10, 10, 10]);
R = diag([10 10]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1 1.4 1.04 1];
exp_param{i}.type = 'lqr'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%31
i = 31;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/2).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([10, 10, 10, 10, 1e14, 1e14]);
R = diag([10 10]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1 1.1 1.01 1];
exp_param{i}.type = 'lqri'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%32
i = 32;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/2).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([10, 10, 10, 10, 1e14, 1e14]);
R = diag([10 10]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1 1.2 1.02 1];
exp_param{i}.type = 'lqri'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%32
i = 32;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/2).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([10, 10, 10, 10, 1e14, 1e14]);
R = diag([10 10]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1 0.9 0.99 1];
exp_param{i}.type = 'lqri'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%33
i = 33;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/2).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([10, 10, 10, 10, 1e14, 1e14]);
R = diag([1000 1000]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1 1.1 1.01 1];
exp_param{i}.type = 'lqri'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%34
i = 34;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/2).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([10, 10, 10, 10, 1e14, 1e14]);
R = diag([1000 1000]);
w = randn(T, 4) .* [.0 .0 .0 .0];
z = randn(T, 2) .* [.0 .0];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1 0.8 0.98 1];
exp_param{i}.type = 'lqri'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%35
i = 35;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/4).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)].*[0.9;0.99] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([10, 10, 10, 10, 1e14, 1e14]);
R = diag([10 10]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];   
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe;
exp_param{i}.type = 'lqgi'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%36 
i = 36;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/4).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)].*[0.9;0.99] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([1, 1, 1, 1, 1e14, 1e14]);
R = diag([1 1]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe;
exp_param{i}.type = 'lqgi'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%37
i = 37;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/4).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)].*[0.9;0.99] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([1, 1, 1, 1, 1e15, 1e15]);
R = diag([1 1]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe;
exp_param{i}.type = 'lqgi'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%38
i = 38;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/4).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)].*[0.9;0.99] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([100, 100, 100, 100, 1e14, 1e14]);
R = diag([100 100]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe;
exp_param{i}.type = 'lqgi'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%39 
i = 39;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/2).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([1, 1, 1, 1, 1e14, 1e14]);
R = diag([1 1]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe;
exp_param{i}.type = 'lqgi'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%40 
i = 40;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/2).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([100, 100, 100, 100, 1e14, 1e14]);
R = diag([100 100]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe;
exp_param{i}.type = 'lqgi'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%41 
i = 41;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/2).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([1, 1, 1, 1]);
R = diag([1 1]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe;
exp_param{i}.type = 'lqg'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%42 
i = 42;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/2).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([100, 100, 100, 100]);
R = diag([100 100]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe;
exp_param{i}.type = 'lqg'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%43
i = 43;
t = (0:0.1:23.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([10, 25, 25, 10]);
R = diag([10 10]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1 0.9 0.99 1];
exp_param{i}.type = 'lqg'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%44
i = 44;
t = (0:0.1:23.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([10, 25, 25, 10]);
R = diag([10 10]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [0.8 0.8 0.98 0.98];
exp_param{i}.type = 'lqg'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%45
i = 45;
t = (0:0.1:23.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([1, 1, 1, 1]);
R = diag([1 1]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [0.8 0.8 0.98 0.98];
exp_param{i}.type = 'lqg'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%46
i = 46;
t = (0:0.1:23.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([100, 100, 100, 100]);
R = diag([100 100]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [0.8 0.8 0.98 0.98];
exp_param{i}.type = 'lqg'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%47
i = 47;
t = (0:0.1:23.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([1000, 1000, 100, 100]);
R = diag([100 100]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [0.8 0.8 0.98 0.98];
exp_param{i}.type = 'lqg'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%48
i = 48;
t = (0:0.1:23.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([1000, 1000, 100, 100]);
R = diag([100 100]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1.2 1.2 1.02 1.02];
exp_param{i}.type = 'lqg'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%49
i = 49;
t = (0:0.1:23.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([1, 1, 1, 1]);
R = diag([1 1]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1.2 1.2 1.02 1.02];
exp_param{i}.type = 'lqg'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%50
i = 50;
t = (0:0.1:23.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([10, 10, 10, 10]);
R = diag([10 10]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1.2 1.2 1.02 1.02];
exp_param{i}.type = 'lqg'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%51 
i = 51;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/2).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([1, 1, 1, 1]);
R = diag([1 1]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe;
exp_param{i}.type = 'lqg'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%52
i = 52;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/2).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([100, 100, 100, 100]);
R = diag([100 100]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe;
exp_param{i}.type = 'lqg'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%53
i = 53;
t = (0:0.1:4.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([10, 25, 25, 10]);
R = diag([10 10]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1 0.9 0.99 1];
exp_param{i}.type = 'lqg'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%54
i = 54;
t = (0:0.1:4.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([10, 25, 25, 10]);
R = diag([10 10]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [0.8 0.8 0.98 0.98];
exp_param{i}.type = 'lqg'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%55
i = 55;
t = (0:0.1:4.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([1, 1, 1, 1]);
R = diag([1 1]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [0.8 0.8 0.98 0.98];
exp_param{i}.type = 'lqg'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%56
i = 56;
t = (0:0.1:4.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([100, 100, 100, 100]);
R = diag([100 100]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [0.8 0.8 0.98 0.98];
exp_param{i}.type = 'lqg'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%57
i = 57;
t = (0:0.1:4.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([1000, 1000, 100, 100]);
R = diag([100 100]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [0.8 0.8 0.98 0.98];
exp_param{i}.type = 'lqg'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%58
i = 58;
t = (0:0.1:4.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([1000, 1000, 100, 100]);
R = diag([100 100]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1.2 1.2 1.02 1.02];
exp_param{i}.type = 'lqg'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%59
i = 59;
t = (0:0.1:4.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([1, 1, 1, 1]);
R = diag([1 1]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1.2 1.2 1.02 1.02];
exp_param{i}.type = 'lqg'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%60
i = 60;
t = (0:0.1:4.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([10, 10, 10, 10]);
R = diag([10 10]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1.2 1.2 1.02 1.02];
exp_param{i}.type = 'lqg'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%61
i = 61;
t = (0:0.1:4.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([10, 10, 10, 10]);
R = diag([10 10]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1.4 1.4 1.04 1.04];
exp_param{i}.type = 'lqg'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%62
i = 62;
t = (0:0.1:4.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([10, 10, 10, 10]);
R = diag([10 10]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [0.6 0.6 0.96 0.96];
exp_param{i}.type = 'lqg'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%63
i = 63;
t = (0:0.1:4.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([10, 10, 10, 10]);
R = diag([10 10]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1 0.6 0.96 1];
exp_param{i}.type = 'lqg'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%64
i = 64;
t = (0:0.1:4.9)'; T = numel(t);
r = ones(2,T).*[xe(2); xe(3)];
Q = diag([10, 10, 10, 10]);
R = diag([10 10]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1 1.4 1.04 1];
exp_param{i}.type = 'lqg'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%65
i = 65;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/2).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([10, 10, 10, 10, 1e14, 1e14]);
R = diag([10 10]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1 1.1 1.01 1];
exp_param{i}.type = 'lqgi'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%66
i = 66;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/2).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([10, 10, 10, 10, 1e14, 1e14]);
R = diag([10 10]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1 1.2 1.02 1];
exp_param{i}.type = 'lqgi'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%67
i = 67;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/2).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([10, 10, 10, 10, 1e14, 1e14]);
R = diag([10 10]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1 0.9 0.99 1];
exp_param{i}.type = 'lqgi'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%68
i = 68;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/2).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([10, 10, 10, 10, 1e14, 1e14]);
R = diag([1000 1000]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1 1.1 1.01 1];
exp_param{i}.type = 'lqgi'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];

%69
i = 69;
t = (0:0.1:23.9)'; T = numel(t);
r = [ones(2,T/4).*[xe(2); xe(3)] ones(2,T/2).*[xe(2); xe(3)].*[1.1;1.01] ones(2,T/4).*[xe(2); xe(3)]];
Q = diag([10, 10, 10, 10, 1e14, 1e14]);
R = diag([1000 1000]);
w = randn(T, 4) .* [.01 .01 .001 .001];
z = randn(T, 2) .* [.01 .001];
exp_param{i}.w = w; exp_param{i}.z = z;
exp_param{i}.t = t; exp_param{i}.r = r; exp_param{i}.x_0 = xe .* [1 0.8 0.98 1];
exp_param{i}.type = 'lqgi'; exp_param{i}.Q = Q; exp_param{i}.R = R; 
exp_param{i}.xout = []; exp_param{i}.yout = []; exp_param{i}.uout = [];