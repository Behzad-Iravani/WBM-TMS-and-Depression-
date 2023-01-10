function [T_boot, CoeffName, idx] = boot_table_logit(T_select,model)


index_select.sham_DEP1 = find(categorical(cellstr(T_select.treatment)) == 'sham'&...
    T_select.Model_Prediction ==0);
index_select.sham_DEP2 = find(categorical(cellstr(T_select.treatment)) == 'sham'&...
    T_select.Model_Prediction ==1);
index_select.active_DEP1 = find(categorical(cellstr(T_select.treatment)) == 'active'&...
    T_select.Model_Prediction ==0);
index_select.active_DEP2 = find(categorical(cellstr(T_select.treatment)) == 'active'&...
    T_select.Model_Prediction ==1);


% ind_boot = @(x) randi(length(x),1*max([length(index_select.sham_DEP1) ...
%                                      length(index_select.sham_DEP2) ...
%                                      length(index_select.active_DEP1) ...
%                                      length(index_select.active_DEP2)]),1);

clear ind_boot
bc = 6; % 15

nboot = round(1e3/(bc*min([length(index_select.sham_DEP1) ...
    length(index_select.sham_DEP2) ...
    length(index_select.active_DEP1) ...
    length(index_select.active_DEP2)])));

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
afterEach(D, @(it) fprintf('%d iteration out of 1000. Elapsed time %1.3f \n', it(1), it(2)))


parfor i = 1:nboot
    
    [T_boot(i,:), CoeffName{i}] = myfitlm([T_select(idx(:,1,i),:)
        T_select(idx(:,2,i),:)
        T_select(idx(:,3,i),:)
        T_select(idx(:,4,i),:)...
        ],model);

send(D,i)

end
end
function [beta_out , name]= myfitlm(T_select,model)
% select out the coefficients that are to be bootstrapped.
flag =1;
try
    beta = fitglm(T_select,model,'Distribution','Binomial', 'Link','logit');
    beta_out  = beta.Coefficients.Estimate;
    if all(abs(beta_out)<30) % converged
        disp('converged')
    else
        disp('NOT converged')
        beta_out = nan(length(beta.Coefficients.Estimate),1);
    end


    name = beta.CoefficientNames;
catch
    beta_out  = nan(1,1);
end
end