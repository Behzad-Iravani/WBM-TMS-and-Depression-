function report_LMM_results2(boot, id)
tmpCI = sign(boot.CI);
index = find(mean(boot.coeff)> 0 & tmpCI(1,:).*tmpCI(2,:) >0);
fprintf(id, 'Significant for DEP1\n')
if ~isempty(index)
    write_line(boot, id , index)
else
    fprintf(id,'none.\n')
end
index = find(mean(boot.coeff)< 0 & tmpCI(1,:).*tmpCI(2,:) >0);
fprintf(id, 'Significant for DEP2\n')
if ~isempty(index)
    write_line(boot, id , index)
else
    fprintf(id,'none.\n')
end
end

function write_line(boot, id , index)
    for i = index
        fprintf(id, boot.name{1}{i});
        fprintf(id, '\nBootStrap: CI = [%1.2f, %1.2f]\n', ...
            boot.CI(1,i), boot.CI(2,i));
    end
end
