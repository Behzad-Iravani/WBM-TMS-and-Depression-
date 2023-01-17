% -*- coding: 'UTF-8' -*-
classdef BHV_LogisticRegression < Results
    % BHV_LogisitcRegression is a subclass of Results. It perform logistic
    % regression the individual items of different measures.
    %   Authors:
    %           Neda Kaboodvand, n.kaboodvand@gmail.com
    %           Behzad Iravani, behzadiravani@gmail.com
    % This function is part of scripts for Macroscopic resting state model predicts
    % theta burst stimulation response: a randomized trial
    % Jan 2023, Stanford, Palo Alto, USA
    %%-------------------------------------------------------------------%%

    properties
        mdl (1,1) struct % structure contains the output of bootstraped logistic regressoion
    end
    methods(Access = public)
        function obj = BHV_LogisticRegression(MonteCarlo, R)
            obj@Results(MonteCarlo);    % call superclass constructor
            obj.TS         = R.TS;
            obj.X          = R.X;
            obj.comorbid   = R.comorbid;
            obj.pdistr     = R.pdistr;
            obj.totalpdistr = R.totalpdistr;
        end

        function obj = run(obj)
            % normalization (zscore) function handel
            nan_zscore = @(x) (x-nanmean(x))./nanstd(x);
            obj.TS.zage = nan_zscore(str2double(cellstr(obj.TS.age)));  % ZSCORE Age
            % add detailed items to the table
            obj.TS = obj.add_details(obj.TS);

            iparm = 0 ;
            number_items = [9, 13, 5];
            % Zscoring ------
            for parms = ["madrs", "cains", "bprs"]
                iparm = iparm +1;
                for items = 1:number_items(iparm)
                    obj.TS.(strcat(parms,sprintf("%d_B", items))) = nan_zscore(...
                        obj.TS.(strcat(parms,sprintf("%d_B", items))));
                    obj.TS.(strcat(parms,sprintf("%d_F", items))) = nan_zscore(...
                        obj.TS.(strcat(parms,sprintf("%d_F", items))));
                end
            end % for parms
            % -----------------
            % BASELINE MODEL
            disp('Baseline ....')
            iparm = 0;
            for parms = ["madrs", "cains", "bprs"]
                iparm = iparm +1;
                model = "Model_Prediction ~ 1 + zage + sex";
                for items = 1:number_items(iparm)
                    model = strcat(model, "+", parms,sprintf("%d_B", items));
                end % end for items
                warning off
                [obj.mdl.baseline.(parms).coeff, obj.mdl.baseline.(parms).name] =...
                    utility.boot_table_logit(obj.TS, model);
                obj.mdl.baseline.(parms).eIDX = find(cellfun(@(x) contains(x,parms), obj.mdl.baseline.(parms).name{1}));
                obj.mdl.baseline.(parms).CI  = quantile(obj.mdl.baseline.(parms).coeff,[.025,.975]);
                disp('plotting!')
                clf
                obj.plot_ColoredColumns(obj.mdl.baseline.(parms))
                title(strcat("\rm", parms))
                print(strcat("result\baseline_", parms,".svg"), '-dsvg')
                % write stats to txt file 
                id = fopen(strcat("result\baseline_", parms, ".txt"),'w');
                obj.report(obj.mdl.baseline.(parms), id);
                fclose(id);
            end % end for parms
            disp('done!')
            %%***********************************************%%
            disp('Follow-up...')
            iparm = 0;
            for parms = ["madrs", "cains", "bprs"]
                iparm = iparm +1;
                model = "Model_Prediction ~ 1 + zage + sex";
                for items = 1:number_items(iparm)
                    model = strcat(model, "+", parms,sprintf("%d_F", items));
                end % end for items
                warning off
                [obj.mdl.followup.(parms).coeff, obj.mdl.followup.(parms).name] =...
                    utility.boot_table_logit(obj.TS, model);
                obj.mdl.followup.(parms).eIDX = find(cellfun(@(x) contains(x,parms), obj.mdl.followup.(parms).name{1}));
                obj.mdl.followup.(parms).CI  = quantile(obj.mdl.followup.(parms).coeff,[.025,.975]);
                disp('plotting!')
                clf
                obj.plot_ColoredColumns(obj.mdl.followup.(parms))
                title(strcat("\rm", parms))
                print(strcat("result\followup_", parms,".svg"), '-dsvg')
                % write stats to txt file 
                id = fopen(strcat("result\followup_", parms, ".txt"),'w');
                obj.report(obj.mdl.followup.(parms), id);
                fclose(id);
            end % end for parms



        end % end run
    end % end public methods

    methods(Access = private)
        function Td = add_details(~, T)
            details = readtable('items.csv');
            Td = [T, details(:,8:end)];
        end % add_details

        function plot_ColoredColumns(~, mdl)
            index = ~cellfun(@(x) contains('zage',x) | contains('sex_F',x) | contains('sex_M',x) | contains('(Intercept)',x), mdl.name{1});
            sig  = sign(mdl.CI(:,index));
            sig  = sig(1,:).*sig(2,:);
            values = nanmean(mdl.coeff(:,index));
            col = [[0 77 129;237 137 93]./255;ones(1,3)./2];
            col =col([2 3 1],:);
            hold on
            b = zeros(1,length(values));
            for i = 1:length(values)
                if values(i) >0 && sig(i)>0
                    b(i) = 1;
                elseif values(i)<0  && sig(i)>0
                    b(i) = -1;
                end
            end
            imagesc(b)
            c = 0;
            for i =find(index)
                c = c + 1;
                if b(c) ==1
                    text(c,1,sprintf('CI = [%1.2f %1.2f]', mdl.CI(1,i),mdl.CI(2,i)),...
                        HorizontalAlignment="center", Rotation= 90 , Color='w')
                else
                    text(c,1,sprintf('CI = [%1.2f %1.2f]', mdl.CI(1,i),mdl.CI(2,i)),...
                        HorizontalAlignment="center", Rotation= 90 )
                end
            end
            axis tight
            line([[1.5:(sum(index)+.5)]'*ones(1,2)]', repmat(ylim(),sum(index),1)' , 'color','k')
            set(gca, 'ytick','', 'xtick', 1:sum(index), 'FontName', 'Arial', 'FontSize', 14)
            colormap(col)
            caxis([-1 1])
            pbaspect([1,.5,1])
        end

        function report(~, boot, id)
            % report writes the bootstrap CI for the logistic regression to
            % text file

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
        

        function write_line(boot, id , index)
            for i = index
                fprintf(id, boot.name{1}{i});
                fprintf(id, '\nBootStrap: CI = [%1.2f, %1.2f]\n', ...
                    boot.CI(1,i), boot.CI(2,i));
            end
        end % write_line 
        end % report 
    end % private methods
end