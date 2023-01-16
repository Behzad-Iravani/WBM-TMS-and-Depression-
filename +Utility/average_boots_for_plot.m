% -*- coding: 'UTF-8' -*-
function T_p = average_boots_for_plot(mdl,T_select,var1, var2)
% average_boots_for_plot avearge over bootstrap iteration for ploting
%   Authors:
%           Neda Kaboodvand, n.kaboodvand@gmail.com
%           Behzad Iravani, behzadiravani@gmail.com
% This function is part of scripts for Macroscopic resting state model predicts
% theta burst stimulation response: a randomized trial
% Jan 2023, Stanford, Palo Alto, USA

%%-----------------------------------------------------------------------%%

T_p = table();
for i = 1:4:size(mdl.(var1).idx,3)
    % sham DEP1
    T_p.group(i)    = categorical({'depression'});
    T_p.treatment(i)= categorical({'sham'});
    T_p.Model_prediction(i) = categorical({'DEP1'});
    T_p.(var2)(i)    = nanmean(T_select.(var2)(mdl.(var1).idx(:,1,i),:));
    % sham DEP2
    T_p.group(i+1)    = categorical({'depression'});
    T_p.treatment(i+1)= categorical({'sham'});
    T_p.Model_prediction(i+1) = categorical({'DEP2'});
    T_p.(var2)(i+1)    = nanmean(T_select.(var2)(mdl.(var1).idx(:,2,i),:));
    % active DEP1
    T_p.group(i+2)    = categorical({'depression'});
    T_p.treatment(i+2)= categorical({'active'});
    T_p.Model_prediction(i+2) = categorical({'DEP1'});
    T_p.(var2)(i+2)    = nanmean(T_select.(var2)(mdl.(var1).idx(:,3,i),:));

    % active DEP2
    T_p.group(i+3)    = categorical({'depression'});
    T_p.treatment(i+3)= categorical({'active'});
    T_p.Model_prediction(i+3) = categorical({'DEP2'});
    T_p.(var2)(i+3)    = nanmean(T_select.(var2)(mdl.(var1).idx(:,4,i),:));

end