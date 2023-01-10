function [T_boot, CoeffName, idx] = boot_table(T_select, index_select,model)

T_select.Model_Prediction = categorical(T_select.Model_Prediction);
index_select.sham_DEP1 = find(categorical(cellstr(T_select.treatment)) == 'sham'&...
    T_select.Model_Prediction =='DEP1');
index_select.sham_DEP2 = find(categorical(cellstr(T_select.treatment)) == 'sham'&...
    T_select.Model_Prediction == 'DEP2');
index_select.active_DEP1 = find(categorical(cellstr(T_select.treatment)) == 'active'&...
    T_select.Model_Prediction =='DEP1');
index_select.active_DEP2 = find(categorical(cellstr(T_select.treatment)) == 'active'&...
    T_select.Model_Prediction == 'DEP2');

nboot = 1e3;
% ind_boot = @(x) randi(length(x),1*max([length(index_select.sham_DEP1) ...
%                                      length(index_select.sham_DEP2) ...
%                                      length(index_select.active_DEP1) ...
%                                      length(index_select.active_DEP2)]),1);

clear ind_boot
bc = 10;
rng(1)
ind_boot(:,:,1) = randi(length(index_select.sham_DEP1),bc*min([length(index_select.sham_DEP1) ...
    length(index_select.sham_DEP2) ...
    length(index_select.active_DEP1) ...
    length(index_select.active_DEP2)]), nboot);
rng(2)
ind_boot(:,:,2) = randi(length(index_select.sham_DEP2),bc*min([length(index_select.sham_DEP1) ...
    length(index_select.sham_DEP2) ...
    length(index_select.active_DEP1) ...
    length(index_select.active_DEP2)]), nboot);
rng(3)
ind_boot(:,:,3) = randi(length(index_select.active_DEP1),bc*min([length(index_select.sham_DEP1) ...
    length(index_select.sham_DEP2) ...
    length(index_select.active_DEP1) ...
    length(index_select.active_DEP2)]), nboot);
rng(4)
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
afterEach(D, @(it) fprintf('%d iteration out of 5000. Elapsed time %1.3f \n', it(1), it(2)))


parfor i = 1:nboot
    % idx(:,1,i) = index_select.sham_DEP1(ind_boot(index_select.sham_DEP1));
    % idx(:,2,i) = index_select.sham_DEP2(ind_boot(index_select.sham_DEP2));
    % idx(:,3,i) = index_select.active_DEP1(ind_boot(index_select.active_DEP1));
    % idx(:,4,i) = index_select.active_DEP2(ind_boot(index_select.active_DEP2));
     t = tic();

    [T_boot(i,:), CoeffName{i}] = myfitlm([T_select(idx(:,1,i),:)
        T_select(idx(:,2,i),:)
        T_select(idx(:,3,i),:)
        T_select(idx(:,4,i),:)...
        ],model);


send(D, [i, toc(t)])
end
end
function [beta_out , name]= myfitlm(T_select,model)
% persistent call_
% if isempty(call_)
%     call_ = 0;
% end
% call_ = call_ +1;
% 
% select out the coefficients that are to be bootstrapped.

try
    beta = fitlm(T_select, model);
    % %     null = regexprep(model, '*', '+');
    % %     betaN = fitlm(T_select, null);


    beta_out  = beta.Coefficients.Estimate;
    name = beta.CoefficientNames;
catch
    beta_out  = nan(1,1);
end
end