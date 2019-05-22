%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamical Modelling and Intelligent Control of Chemical Reactive Systems
% 
% Author: Otacilio Bezerra Leite Neto
%
%   exo_dynamics.m
%       This script contains the instructions for running and visualize experiments 
%       over the dynamics of the Exothermal Continuous-Stirred Tank (CSTR) system.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preamble %%
clc; close all; clear all;
% cd /home/minho/Documents/Minho/Control_ChemPlants/
% cd /home/minhotmog/Dropbox/Research/TCC/Codes/

% Sets the Default Rendere to tbe the Painters
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
        91,192,222    % 7 - Bootstrap Light Blue
        92,184,92     % 8 - Bootstrap Green
        66,139,202    % 9 - Bootstrap Blue
        255,167,0     % 10 - Google Yellow    
       ]/255;  
                          
%% %%%%%%%%%%%%
%  EXOTHERMAL CSTR %%
%  %%%%%%%%%%%%%%
%% Model Loading %%
run exothermal_cstr/exo_model.m
%load('data/exo_cstr_model.mat')
