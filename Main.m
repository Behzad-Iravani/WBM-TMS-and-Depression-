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

% PLOTING CONNECTOME 
con = cellfun(@simple_comp_FC,Empirical.TC,'UniformOutput',false);
FC1 = mean(cat(3,con{Results.TS.Model_Prediction == 1}),3).*~eye(68); % 
FC2 = mean(cat(3,con{Results.TS.Model_Prediction == 2}),3).*~eye(68); %
subplot(121)
imagesc(FC1)
colorbar()
axis square
axis off
axis tight
title('DEP1 -- empirical FC')
subplot(122)
imagesc(FC2)
colorbar()
axis square
axis off
axis tight
title('DEP2 -- empirical FC')

% % % ploting simulated FCs for 
% % subplot(223)
% % imagesc(Results.X.dep1.FC)
% % axis square
% % axis off
% % title('DEP1 -- simulated FC')
% % subplot(224)
% % imagesc(Results.X.dep2.FC)
% % axis square
% % axis off
% % title('DEP2 -- simulate FC')
% % 
% fprintf('>> SDEP1 and SDEP2 FC similarity = %1.2f\n',...
%     corr(squareform(Results.X.dep1.FC)', squareform(Results.X.dep2.FC)'))
fprintf('>> EDEP1 and EDEP2 FC similarity = %1.2f\n',...
    corr(squareform(FC1)', squareform(FC2)'))
print result\FC_dep1_dep2.png -r300 -dpng
% fprintf('>> SDEP1 and EDEP1 FC similarity = %1.2f\n',...
%     corr(squareform(Results.X.dep1.FC)', squareform(FC1)'))
% fprintf('>> SDEP2 and EDEP2 FC similarity = %1.2f\n',...
%     corr(squareform(Results.X.dep2.FC)', squareform(FC2)'))
% fprintf('>> SDEP1 and EDEP2 FC similarity = %1.2f\n',...
%     corr(squareform(Results.X.dep1.FC)', squareform(FC2)'))
% fprintf('>> SDEP2 and EDEP1 FC similarity = %1.2f\n',...
%     corr(squareform(Results.X.dep2.FC)', squareform(FC1)'))


con = cellfun(@simple_comp_FC,Empirical.TC,'UniformOutput',false);
FC1Sham = mean(cat(3,con{Results.TS.Model_Prediction == 1 & categorical(cellstr(Results.TS.treatment)) == 'sham'}),3).*~eye(68); % 2 more populated is DEP1
FC2Sham = mean(cat(3,con{Results.TS.Model_Prediction == 2& categorical(cellstr(Results.TS.treatment)) == 'sham'}),3).*~eye(68); % 1 less populated is DEP2
FC1Active = mean(cat(3,con{Results.TS.Model_Prediction == 1& categorical(cellstr(Results.TS.treatment)) == 'active'}),3).*~eye(68); % 1 less populated is DEP2
FC2Active = mean(cat(3,con{Results.TS.Model_Prediction == 2& categorical(cellstr(Results.TS.treatment)) == 'active'}),3).*~eye(68); % 1 less populated is DEP2

subplot(221)
imagesc(FC1Sham)
colorbar()
clim([0,1.2])
axis square
axis off
axis tight                                                                  
title('DEP1 Sham -- empirical FC', ['n = ', num2str(sum(Results.TS.Model_Prediction == 1 & categorical(cellstr(Results.TS.treatment)) == 'sham'))])
subplot(222)
imagesc(FC2Sham)
colorbar()
clim([0,1.2])
axis square
axis off
axis tight
title('DEP2 Sham -- empirical FC', ['n = ', num2str(sum(Results.TS.Model_Prediction == 2 & categorical(cellstr(Results.TS.treatment)) == 'sham'))])
subplot(223)
imagesc(FC1Active)
colorbar()
clim([0,1.2])
axis square
axis off
axis tight                                                                   
title('DEP1 Active -- empirical FC', ['n = ', num2str(sum(Results.TS.Model_Prediction == 1 & categorical(cellstr(Results.TS.treatment)) == 'active'))]) 
subplot(224)
imagesc(FC2Active)
colorbar()
clim([0,1.2])
axis square
axis off
axis tight
title('DEP2 Active -- empirical FC',['n = ', num2str(sum(Results.TS.Model_Prediction == 2 & categorical(cellstr(Results.TS.treatment)) == 'active'))]) 
print result\FC_dep1_dep2&treatments.png -r300 -dpng


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
%% 3.2	Individual factors of MADRS-S, BPRS AFF and CAINS differentiate depression subtypes at baseline 
LG = BHV_LogisticRegression({MonteCarloJob_30, MonteCarloJob_50, MonteCarloJob_70}, Results);
% run bootstrapped linear regression
LG = LG.run();
save LG.mat LG

%% 3.3	The therapeutic effect of iTBS on dmPFC is modulated by intrinsic brain dynamics at baseline 
LM = BHV_LinearRegression({MonteCarloJob_30, MonteCarloJob_50, MonteCarloJob_70}, Results);
% run bootstrapped linear regression
LM = LM.run();
% ploting 
warning off
LM.plot
warning on
print -dpng -r300 result\ViolinMeasure.png
print -dsvg -vector result\ViolinMeasure.svg
save LM.mat LM
% % % bar([mean(Results.TS.MetaStability(Results.TS.Model_Prediction == 1))
% % %     mean(Results.TS.MetaStability(Results.TS.Model_Prediction == 2))])
% % %
% % % bar([mean(Results.TS.Synchrony(Results.TS.Model_Prediction == 1))
% % %     mean(Results.TS.Synchrony(Results.TS.Model_Prediction == 2))])

%% 3.4	The iTBS treatment modulates items in a certain depression subgroup 
itANOVA = itemsANOVA({MonteCarloJob_30, MonteCarloJob_50, MonteCarloJob_70}, Results, LG);
itANOVA.run
save itANOVA.mat itANOVA
% Violin/box ploting 
clf
ns = factor(length(fieldnames(itANOVA.stat))-1);
iplot = 0;
for f = fieldnames(itANOVA.stat)'
    if ~strcmp(f{:}, 'id')
        writecell(itANOVA.ANOVAmdl.(f{:}).tbl, strcat('result\', f{:}, '.csv'))


        iplot = iplot +1;
        if length(ns) == 2
            subplot(ns(1), ns(2), iplot)
        else
        subplot(1, ns, iplot)
        end
        itANOVA.plot_test_violin(itANOVA.TS, f{:})
        title(regexprep(strcat('\rm', f{:}), '(delta)', '\\Delta'))
        set(gca,'FontName', 'Arial', 'FontSize', 12)
    end % if
end % for
print -dpng -r300 result\ItemsViolinMeasure.png
print -dsvg result\ItemsViolinMeasure.svg

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





