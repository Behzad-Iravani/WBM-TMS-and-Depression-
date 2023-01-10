% -*- coding: 'UTF-8' -*-
% Monte Carlo Job
clc
clear all
%% load Emorical data
load Data\Empirical.mat
% run and save monte carlo simulations
runMot(Empirical, 20) % The first 500 were run with default seed 15

