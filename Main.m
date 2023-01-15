% -*- coding: 'UTF-8' -*-
%%  Whole brain modeling project of depression 

%{
This is the main script to run the analysis of the manuscript entilied as
"Macroscopic resting state model predicts theta burst stimulation response:
a randomized trial"

Authors:
Neda Kaboodvand 
n.kaboodvand@gmail.com 
Behzad Iravani
behzadiravani@gmail.com


Palo Alto, june 2022
revised jan 2023
%}
%%-----------------------------------------------------------------------%%
clear all, clc

%% loading the function data 
load Data\T_plusBehavior.mat
% find DEP subjects
T = T(~isnan(T.fMRIses1) & categorical(cellstr(T.group)) == 'depression',:);

subjs = dir('TCs_subjDK\ses-baseline\subj*');
c = 0;
for sub = T.ID'
    c = c + 1;
    fprintf('load TCs for %s\n', sub{:})
%     load(strcat("data\",sub{:},filesep, "Atlas_wmparcTCs_ses1.mat"))
    ind = find(cellfun(@(x) contains(x,sub{:}), {subjs.name}));
    load(fullfile(subjs(ind).folder, subjs(ind).name, 'TCs.mat'))
    load(fullfile(subjs(ind).folder, subjs(ind).name, 'ROIs.mat'))
    
    % extract cortical time courses
    TC{c}   = utility.find_cortical_nodes(TCs, ROIs);
    FQ(c,:) = utility.find_freq(TC{c}); 
end
disp('Average the empirical data')
Empirical.treatment = T.treatment;
Empirical.TC    =  TC;
Empirical.Conn  = arrayfun(@(j) corrcoef(Empirical.TC{j}),1:length(Empirical.TC),'UniformOutput',false);
% frequecny of BOLD activity
Empirical.freq  = median(FQ);
% fsample
Empirical.fsamp = 1/2; % 1 over TR 

% Construct the model
Model = WBM;
% Intiate the model parameters
Model.intiate_model
% add empirical data
Model.Empirical = Empirical;
Model.SC        = Model.SC./repmat(sum(Model.SC,2), 1,68);
%% Revison search

Model.search.A = 0;%single([-.2:.01:.2]); 
Model.search.G = 0;%single([0.005:0.005:0.05]);
Model.search.F = 0;%single([-1:0.1:0]);
Model.search.M = 0;%single([0.1:0.03:.5]);
%% Grid search and simulation 
% simulating and saving
% Model.simulate()
%% running Monte Carlo permuatations
MonteCarloJob = MonteCarlo


