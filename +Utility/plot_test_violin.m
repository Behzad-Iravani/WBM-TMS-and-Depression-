function plot_test_violin(T_detailed_selected, field)
dat = [T_detailed_selected.(field)(categorical(cellstr(T_detailed_selected.treatment)) == "active" & ...
    categorical(T_detailed_selected.Model_Prediction) =="DEP1") - ...
    nanmedian(T_detailed_selected.(field)(categorical(cellstr(T_detailed_selected.treatment)) == "sham" & ...
    categorical(T_detailed_selected.Model_Prediction) =="DEP1"))];

V = Utility.Violin(dat,1,'BoxWidth', .05, 'BoxColor', zeros(1,3));
V.BoxPlot.FaceAlpha = .6;
dat = [T_detailed_selected.(field)(categorical(cellstr(T_detailed_selected.treatment)) == "active" & ...
    categorical(T_detailed_selected.Model_Prediction) =="DEP2")-...
    nanmedian(T_detailed_selected.(field)(categorical(cellstr(T_detailed_selected.treatment)) == "sham" & ...
    categorical(T_detailed_selected.Model_Prediction) =="DEP2"))];
col = lines(4);
V = Utility.Violin(dat,2,'BoxWidth', .05, 'BoxColor', zeros(1,3));
V.ViolinPlot.FaceColor = col(4,:);
V.ScatterPlot.MarkerFaceColor = col(4,:);
V.BoxPlot.FaceAlpha = .6;

% ylim(max(abs(ylim)).*[-1,1])
ylim([3].*[-1,1]);
set(gca, 'YTick',[-3,0,3])
set(gca, 'XTick',[1,2,4,5],'XTickLabel',["DEP1", "DEP2"])
axis square
drawnow
end