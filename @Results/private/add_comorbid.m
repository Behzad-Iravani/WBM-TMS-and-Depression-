% -*- coding: 'UTF-8' -*-
function T = add_comorbid(~,path2csv, T)
% add_comorbid is a private method of Results. it adds the comorbid data to
% the table TS for chi-squared analysis and to assess the realtionship
% between model prediction and comorbidity.
%   Authors:
%           Neda Kaboodvand, n.kaboodvand@gmail.com
%           Behzad Iravani, behzadiravani@gmail.com
% This function is part of scripts for Macroscopic resting state model predicts
% theta burst stimulation response: a randomized trial
% Jan 2023, Stanford, Palo Alto, USA 

CoMorbidity     = readtable(path2csv);

[IndexC,IndexT] = find(CoMorbidity.PatientNo_ ==...
    cellfun(@str2num,cellfun(@(x) regexp(x,'\d+', 'match'),T.id))');
for CName  = CoMorbidity.Properties.VariableNames(3:7)
    T.(CName{:})          = nan(height(T),1);
    T.(CName{:})(IndexT)  = CoMorbidity.(CName{:})(IndexC);
end
end