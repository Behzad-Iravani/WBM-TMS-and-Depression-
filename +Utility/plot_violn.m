function plot_violn(T_select,fieldname)
persistent call
    if isempty(call)
        call = 0;
    end
nanz = @(x) (x-nanmean(x))./nanstd(x);
T_select.treatment = categorical(T_select.treatment);
T_select.Model_prediction = categorical(T_select.Model_prediction);

% line([1,2], [nanmedian(T_select.(fieldname)(categorical(cellstr(T_select.treatment)) == "sham" & ...
%     T_select.Model_prediction =="DEP1")),...
%     nanmedian(T_select.(fieldname)(categorical(cellstr(T_select.treatment)) == "sham" & ...
%     T_select.Model_prediction =="DEP2"))]...
%     , 'LineStyle', '-', 'Color','k', 'LineWidth', 1.75)
T_select.(fieldname) = nanz(T_select.(fieldname));
V = utility.Violin(T_select.(fieldname)(T_select.treatment == "sham" & ...
    T_select.Model_prediction =="DEP1"),1,'BoxWidth', .05, 'BoxColor', zeros(1,3),'ShowMean',true);
V.MedianPlot.SizeData = V.MedianPlot.SizeData/2;
V.ScatterPlot.SizeData    = 5;
V.ScatterPlot.MarkerFaceAlpha = .01;
V.ScatterPlot.MarkerEdgeAlpha = .1;
V.ScatterPlot.MarkerEdgeColor = V.ScatterPlot.MarkerFaceColor;

V = utility.Violin(T_select.(fieldname)(T_select.treatment == "sham" & ...
    T_select.Model_prediction =="DEP2"),2,'BoxWidth', .05, 'BoxColor', zeros(1,3),'ShowMean',true);
V.MedianPlot.SizeData = V.MedianPlot.SizeData/2;
V.ScatterPlot.SizeData    = 5;
V.ScatterPlot.MarkerFaceAlpha = .01;
V.ScatterPlot.MarkerEdgeAlpha = .1;
V.ScatterPlot.MarkerEdgeColor = V.ScatterPlot.MarkerFaceColor;

% line([4,5], [nanmedian(T_select.(fieldname)(categorical(cellstr(T_select.treatment)) == "active" & ...
%     T_select.Model_prediction =="DEP1")),...
%     nanmedian(T_select.(fieldname)(categorical(cellstr(T_select.treatment)) == "active" & ...
%     T_select.Model_prediction =="DEP2"))]...
%     , 'LineStyle', '-', 'Color','k', 'LineWidth', 1.75)
V= utility.Violin(T_select.(fieldname)(T_select.treatment == "active" & ...
    T_select.Model_prediction =="DEP1"),4,'BoxWidth', .05, 'BoxColor', zeros(1,3),'ShowMean',true);
V.MedianPlot.SizeData = V.MedianPlot.SizeData/2;
V.ScatterPlot.MarkerFaceAlpha = .01;
V.ScatterPlot.SizeData    = 5;
V.ScatterPlot.MarkerEdgeAlpha = .1;
V.ScatterPlot.MarkerEdgeColor = V.ScatterPlot.MarkerFaceColor;

V = utility.Violin(T_select.(fieldname)(T_select.treatment == "active" & ...
    T_select.Model_prediction =="DEP2"),5,'BoxWidth', .05, 'BoxColor', zeros(1,3),'ShowMean',true);
V.MedianPlot.SizeData = V.MedianPlot.SizeData/2;
V.ScatterPlot.MarkerFaceAlpha = .01;
V.ScatterPlot.SizeData    = 5;
V.ScatterPlot.MarkerEdgeAlpha = .1;
V.ScatterPlot.MarkerEdgeColor = V.ScatterPlot.MarkerFaceColor;
drawnow
% ylim(max(abs(ylim)).*[-1,1])
ylim([3].*[-1,1]);
set(gca, 'YTick',[-3,0,3])
if call == 0
set(gca, 'XTick',[1,2,4,5],'XTickLabel',["Sham:DEP1", "Sham:DEP2" ,"Active:DEP1" ,"Active:DEP2"])

else
    set(gca, 'XTick',[1,2,4,5],'XTickLabel',['' '' '' ''])
end

call = call + 1;
axis square
end