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

Model.search.A = single([-.2:.01:.2]);
Model.search.G = single([0.005:0.005:0.05]);
Model.search.F = single([-1:0.1:0]);
Model.search.M = single([0.1:0.03:.5]);
%% Grid search and simulation
% simulating and saving
% Model.simulate()
%% Running Monte Carlo permuatations
MonteCarloJob_30 = MonteCarlo({fullfile(pwd(), 'Data\simulated.dat')},.... path to .dat simulated file
    .3,... threshold
    [],... numberofperm
    [],... EachItofSim
    []... random seed (default value for reproducibility)
    );
% MonteCarloJob_30.run(Empirical.TC)
%%-------------------------------------------------------------------------%%
MonteCarloJob_50 = MonteCarlo({fullfile(pwd(), 'Data\simulated.dat')},.... path to .dat simulated file
    1-.5,... threshold
    [],... numberofperm
    [],... EachItofSim
    []... random seed (default value for reproducibility)
    );
% MonteCarloJob_50.run(Empirical.TC)
%%-------------------------------------------------------------------------%%
MonteCarloJob_70 = MonteCarlo({fullfile(pwd(), 'Data\simulated.dat')},.... path to .dat simulated file
    .7,... threshold
    [],... numberofperm
    [],... EachItofSim
    []... random seed (default value for reproducibility)
    );
% MonteCarloJob_70.run(Empirical.TC)
%% Results
clear Results
Results = Results({MonteCarloJob_30, MonteCarloJob_50, MonteCarloJob_70});
% plot distributions of parameters for every thershold.
Results = Results.plot_parameters_distributions(Model.search);
% ------------------------------------------------------
Results = Results.plot_synchrony_metastability; % also plots spiderplot for parameters;
close all
% ------------------------------------------------------
Results.plot_agg_dist(); % plot aggergate thresholds (i.e., .3,.5,.7)

% Stratifying patients
Results = Results.Stratifying(Empirical, T);

fprintf('>> N(DEP1) = %d and N(DEP2) = %d\n', sum(Results.TS.Model_Prediction == 1), sum(Results.TS.Model_Prediction == 2))
% chi squared test of model prediction and treatment 
[~,chi2,p_value] = crosstab(Results.TS.Model_Prediction, Results.TS.treatment);
fprintf('chi-squared test model prediction and teatment = %1.2f, p = %1.2f\n',...
   chi2, p_value)
%%********** MOVEMENT AND MOEDL PRED *********************%% 

%%********** COMORBIDITY AND MODEL PRED ******************%%
Results.chi2_comorbid('CAINS_DEP_220113.csv')
for f = fieldnames(Results.comorbid.chi2)'
fprintf('>> %s = %1.2f, p > %1.2f\n', f{:}, Results.comorbid.chi2.(f{:}), Results.comorbid.p.(f{:}))
end
%% 3.2	The therapeutic effect of iTBS on dmPFC is modulated by intrinsic brain dynamics at baseline 
LM = BHV_LinearRegression({MonteCarloJob_30, MonteCarloJob_50, MonteCarloJob_70}, Results);
% run bootstrapped linear regression
LM = LM.run();
% ploting 
LM.plot
print -dpng -r300 result\ViolinMeasure.png
save LM.mat LM
% % % bar([mean(Results.TS.MetaStability(Results.TS.Model_Prediction == 1))
% % %     mean(Results.TS.MetaStability(Results.TS.Model_Prediction == 2))])
% % %
% % % bar([mean(Results.TS.Synchrony(Results.TS.Model_Prediction == 1))
% % %     mean(Results.TS.Synchrony(Results.TS.Model_Prediction == 2))])
%% 3.3	Individual factors of MADRS-S, BPRS AFF and CAINS differentiate depression subtypes at baseline 


%% Plot simulated TCs for major brain lobes
load('+utility\cortical_ROIs.mat')
sp = [1,2,5,6];
i = 0;
clear p
close all
figure('Units','normalized',Position=[0,0,.75,.5])
ETC = cat(3,TC{:});
for brain_lobes = ["Frontal", "Temporal", "Parietal", "Occipital"]
    i = i+1;
    subplot(2,4,sp(i))

    p(1) = plot(zscore(mean(Results.X.dep1.TC(:,strcmp(cortical_ROIs(:,2), brain_lobes)),2)), 'Color',[0 77 129]./255);
    hold on
    p(2) = plot(zscore(mean(Results.X.dep2.TC(:,strcmp(cortical_ROIs(:,2), brain_lobes)),2)), 'Color', [237 137 93]./255);
    hold off
    box off
    if i==1
        xlabel('Time')
        legend(p,{'Sim-DEP1', 'Sim-DEP2'},'Box', 'off')
    end
    title(brain_lobes)
    ylim([-3,3])
end
i = 0;
sp = [3,4,7,8];
clear p
for brain_lobes = ["Frontal", "Temporal", "Parietal", "Occipital"]
    i = i+1;
    subplot(2,4,sp(i))


    p(1) = plot(zscore(squeeze(...
        mean(...
        mean(...
        ETC(:,strcmp(cortical_ROIs(:,2), brain_lobes),Results.TS.Model_Prediction == 1),3)...
        ,2))), 'Color',[0 77 129]./255, 'LineStyle','-');
    hold on
    p(2) = plot(zscore(squeeze(...
        mean(...
        mean(...
        ETC(:,strcmp(cortical_ROIs(:,2), brain_lobes),Results.TS.Model_Prediction == 2),3)...
        ,2))), 'Color',[237 137 93]./255, 'LineStyle','-');

    hold off
    box off
    if i==1
        xlabel('Time')
        legend(p,{'DEP1', 'DEP2'},'Box', 'off')
    end
    title(brain_lobes)
    ylim([-3,3])
end

print result\BrainLobes_simulatedTCs.png -dpng





