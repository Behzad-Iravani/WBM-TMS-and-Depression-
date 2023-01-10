function TC_cortical = find_cortical_nodes(TCs, ROIs_ind)
load("+utility\ROIS.mat")
load("+utility\cortical_ROIS.mat")
% load("data\Atlas_wmparc_Labels.mat")

% find the cortical node indices
[~,j] = find(ROIs == ROIs_ind');
TC_cortical = TCs(:,j);

end