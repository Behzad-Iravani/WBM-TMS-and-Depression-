% -*- coding: 'UTF-8' -*-
function [T_boot, CoeffName, idx] = boot_table(T_select,model)
% boot_table runs a bootstrped linear regession model.
% Inputs:
%           T_select        - Table (data)
%           model           - a string that defines the linear model
%   Authors:
%           Neda Kaboodvand, n.kaboodvand@gmail.com
%           Behzad Iravani, behzadiravani@gmail.com
% This function is part of scripts for Macroscopic resting state model predicts
% theta burst stimulation response: a randomized trial
% Jan 2023, Stanford, Palo Alto, USA

%%------------------------------------------------------------------------%%

T_select.Model_Prediction = categorical(T_select.Model_Prediction);% converts to categorical variable 
% extracting the indcies for treatment and dep subtypes
index_select.sham_DEP1 = find(categorical(cellstr(T_select.treatment)) == 'sham'&...
    T_select.Model_Prediction =='DEP1');
index_select.sham_DEP2 = find(categorical(cellstr(T_select.treatment)) == 'sham'&...
    T_select.Model_Prediction == 'DEP2');
index_select.active_DEP1 = find(categorical(cellstr(T_select.treatment)) == 'active'&...
    T_select.Model_Prediction =='DEP1');
index_select.active_DEP2 = find(categorical(cellstr(T_select.treatment)) == 'active'&...
    T_select.Model_Prediction == 'DEP2');
% ----
nboot = 1e3; % number of boostraping 
clear ind_boot
bc = 10;     % defines the minimum number of resampling with replacement for each iterations 
rng(1)       % setting the random seed for reporducibility  
ind_boot(:,:,1) = randi(length(index_select.sham_DEP1),bc*min([length(index_select.sham_DEP1) ...
    length(index_select.sham_DEP2) ...
    length(index_select.active_DEP1) ...
    length(index_select.active_DEP2)]), nboot);
rng(2)       % setting the random seed for reporducibility
ind_boot(:,:,2) = randi(length(index_select.sham_DEP2),bc*min([length(index_select.sham_DEP1) ...
    length(index_select.sham_DEP2) ...
    length(index_select.active_DEP1) ...
    length(index_select.active_DEP2)]), nboot);
rng(3)      % setting the random seed for reporducibility
ind_boot(:,:,3) = randi(length(index_select.active_DEP1),bc*min([length(index_select.sham_DEP1) ...
    length(index_select.sham_DEP2) ...
    length(index_select.active_DEP1) ...
    length(index_select.active_DEP2)]), nboot);
rng(4)      % setting the random seed for reporducibility
ind_boot(:,:,4) = randi(length(index_select.active_DEP2),bc*min([length(index_select.sham_DEP1) ...
    length(index_select.sham_DEP2) ...
    length(index_select.active_DEP1) ...
    length(index_select.active_DEP2)]), nboot);

idx = permute(ind_boot,[1,3,2]);
idx = cat(2,index_select.sham_DEP1(idx(:,1,:)),...
    index_select.sham_DEP2(idx(:,2,:)),...
    index_select.active_DEP1(idx(:,3,:)),...
    index_select.active_DEP2(idx(:,4,:)));
clear myfitlm
D = parallel.pool.DataQueue;
afterEach(D, @(it) fprintf('%d iteration out of 1000. Elapsed time %1.3f \n', it(1), it(2)))

parfor i = 1:nboot

    t = tic();

    [T_boot(i,:), CoeffName{i}] = myfitlm([T_select(idx(:,1,i),:)
        T_select(idx(:,2,i),:)
        T_select(idx(:,3,i),:)
        T_select(idx(:,4,i),:)...
        ],model);


    send(D, [i, toc(t)])
end
end
function [beta_out , name]= myfitlm(T_select,model) % LME function
try
    beta = fitlm(T_select, model);
    beta_out  = beta.Coefficients.Estimate;
    name = beta.CoefficientNames;
catch
    beta_out  = nan(1,1);
end
end