function stat = ttest2_phenotype(T_select,fieldname)

[~, stat.sham.p, stat.sham.CI, stat.sham.stat] = ttest2(T_select.(fieldname)(categorical(cellstr(T_select.treatment)) == "sham" & ...
    categorical(T_select.Model_Prediction) =="DEP1"),...
T_select.(fieldname)(categorical(cellstr(T_select.treatment)) == "sham" & ...
    categorical(T_select.Model_Prediction) =="DEP2"));


[~, stat.active.p, stat.active.CI, stat.active.stat] = ttest2(T_select.(fieldname)(categorical(cellstr(T_select.treatment)) == "active" & ...
    categorical(T_select.Model_Prediction) =="DEP1"),...
T_select.(fieldname)(categorical(cellstr(T_select.treatment)) == "active" & ...
    categorical(T_select.Model_Prediction) =="DEP2"));



X1 = T_select.(fieldname)(categorical(cellstr(T_select.treatment)) == "active" & ...
    categorical(T_select.Model_Prediction) =="DEP1") - ...
  nanmedian(T_select.(fieldname)(categorical(cellstr(T_select.treatment)) == "sham" & ...
    categorical(T_select.Model_Prediction) =="DEP1"));
X2 = T_select.(fieldname)(categorical(cellstr(T_select.treatment)) == "active" & ...
    categorical(T_select.Model_Prediction) =="DEP2")-...
      nanmedian(T_select.(fieldname)(categorical(cellstr(T_select.treatment)) == "sham" & ...
    categorical(T_select.Model_Prediction) =="DEP2"));


% [~, stat.multi.p, stat.multi.CI, stat.multi.stat] = ttest2(X1 , X2...
% , 'Vartype','unequal');

[stat.multi.p, stat.multi.CI, stat.multi.stat] = bootest(X1 , X2);

end
function [p, CI, stat] = bootest(x,y)
rng(1)
index1 = randi(length(x),round(3*length(x)),3e2);
rng(2)
index2 = randi(length(y),round(3*length(y)),3e2);
for ii = 1:3e2
      [~,tmpP(ii), tmpCI(:,ii), tmpStat(ii)]= ttest2(x(index1(:,ii)), y(index2(:,ii)));
end
p = mean(tmpP);
CI = mean(tmpCI,2);
Xt = [tmpStat.tstat];
stat.tstat = mean(Xt(~isinf(Xt)));
stat.df    = mean([tmpStat.df]);

end


