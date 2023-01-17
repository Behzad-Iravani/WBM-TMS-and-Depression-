% -*- coding: 'UTF-8' -*-
classdef itemsANOVA < BHV_LogisticRegression
    %  itemsANOVA is a subclass of Results and explores the effect of treatment
    %  in cretain depression subgroup using hierarchical ANOVA
    %   Authors:
    %           Neda Kaboodvand, n.kaboodvand@gmail.com
    %           Behzad Iravani, behzadiravani@gmail.com
    % This function is part of scripts for Macroscopic resting state model predicts
    % theta burst stimulation response: a randomized trial
    % Jan 2023, Stanford, Palo Alto, USA

    properties
        ANOVAmdl (1,1) struct % structure contains the ouput of the nesting ANOVA
        stat (1,1) struct % structure contains the ttest statistic
    end

    methods(Access = public)
        function obj = itemsANOVA(MonteCarlo, R, LG) % constructor method
            obj@BHV_LogisticRegression(MonteCarlo, R)
            obj.TS = LG.TS;
        end % end constructor

        function run(obj)
             M = containers.Map([1,2],{'DEP1', 'DEP2'});
            obj.TS.Model_Prediction = arrayfun(@(j) M(j),obj.TS.Model_Prediction,'UniformOutput',0);

            iparm = 0 ; % parms counter
            number_items = [9, 13, 5];
            icounter = 0;   % parms and items counter
            sigCounter = 0; % significant counter
            for parms = ["madrs", "cains", "bprs"]
                iparm = iparm +1;
                for items = 1:number_items(iparm)
                    icounter = icounter + 1;
                    obj.TS.(strcat("delta",parms, string(num2str(items)))) = obj.computedelat(obj.TS, strcat(parms, string(num2str(items)))); % computed delta

                    [~,obj.ANOVAmdl.(strcat("delta",parms, string(num2str(items)))).tbl{items},...
                        obj.ANOVAmdl.(strcat("delta",parms, string(num2str(items)))).ANOVA{items}]...
                        = anovan(obj.TS.(strcat("delta",parms, string(num2str(items)))),...
                        {obj.TS.zage,...
                        obj.TS.sex,...
                        obj.TS.treatment,...
                        obj.TS.Model_Prediction}, "model",[1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1],'nested',[0,0,0,0
                        0 0 0 0
                        0 0 0 1
                        0 0 0 0], "varnames",{'age','sex','treatment','DEP'});

                    obj.ANOVAmdl.(strcat("delta",parms, string(num2str(items)))).p = obj.ANOVAmdl.(strcat("delta",parms, string(num2str(items)))).tbl{items}{4,7};
                    obj.ANOVAmdl.(strcat("delta",parms, string(num2str(items)))).F = obj.ANOVAmdl.(strcat("delta",parms, string(num2str(items)))).tbl{items}{4,6};
                    
                    if obj.ANOVAmdl.(strcat("delta",parms, string(num2str(items)))).p<=.065
                        sigCounter = sigCounter + 1;
                        obj.stat.id(sigCounter) = icounter;
                        obj = obj.ttest2_phenotype(strcat("delta",parms, string(num2str(items))));
                    end %if
                end % items
            end % parms
            close all hidden
        end % run

        function obj = ttest2_phenotype(obj,fieldname)
           
            [~, obj.stat.(fieldname).sham.p, obj.stat.(fieldname).sham.CI, obj.stat.(fieldname).sham.stat] = ttest2(obj.TS.(fieldname)(categorical(cellstr(obj.TS.treatment)) == "sham" & ...
                categorical(obj.TS.Model_Prediction) =="DEP1"),...
                obj.TS.(fieldname)(categorical(cellstr(obj.TS.treatment)) == "sham" & ...
                categorical(obj.TS.Model_Prediction) =="DEP2"));


            [~, obj.stat.(fieldname).active.p, obj.stat.(fieldname).active.CI, obj.stat.(fieldname).active.stat] = ttest2(obj.TS.(fieldname)(categorical(cellstr(obj.TS.treatment)) == "active" & ...
                categorical(obj.TS.Model_Prediction) =="DEP1"),...
                obj.TS.(fieldname)(categorical(cellstr(obj.TS.treatment)) == "active" & ...
                categorical(obj.TS.Model_Prediction) =="DEP2"));



            X1 = obj.TS.(fieldname)(categorical(cellstr(obj.TS.treatment)) == "active" & ...
                categorical(obj.TS.Model_Prediction) =="DEP1") - ...
                nanmedian(obj.TS.(fieldname)(categorical(cellstr(obj.TS.treatment)) == "sham" & ...
                categorical(obj.TS.Model_Prediction) =="DEP1"));
            X2 = obj.TS.(fieldname)(categorical(cellstr(obj.TS.treatment)) == "active" & ...
                categorical(obj.TS.Model_Prediction) =="DEP2")-...
                nanmedian(obj.TS.(fieldname)(categorical(cellstr(obj.TS.treatment)) == "sham" & ...
                categorical(obj.TS.Model_Prediction) =="DEP2"));

            [obj.stat.(fieldname).multi.p, obj.stat.(fieldname).multi.CI, obj.stat.(fieldname).multi.stat] = bootest(X1 , X2);

        
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
        end % bootest 
        end % ttest2_phenotype
    end % public methods
    methods(Static, Access = private)
        function delta = computedelat(Table, parameter)
            delta = Table.(strcat(parameter,"_F"))- Table.(strcat(parameter,"_B"));
        end

    end % Static and private methods

end