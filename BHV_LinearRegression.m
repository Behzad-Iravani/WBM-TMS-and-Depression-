% -*- coding: 'UTF-8' -*-
classdef BHV_LinearRegression < Results
    % BHV_LinearRegression is a subclass of Results. It performs linear
    % regression on the behavioral data and assesses the association
    % between the identfied sub-cohorts and behavior
    %   Authors:
    %           Neda Kaboodvand, n.kaboodvand@gmail.com
    %           Behzad Iravani, behzadiravani@gmail.com
    % This function is part of scripts for Macroscopic resting state model predicts
    % theta burst stimulation response: a randomized trial
    % Jan 2023, Stanford, Palo Alto, USA

    properties
        mdl (1,1) struct % structure contains the output of bootstraped linear regressoion
    end
    methods(Access = public)
        function obj = BHV_LinearRegression(MonteCarlo, R)
            obj@Results(MonteCarlo);    % call superclass constructor
            obj.TS         = R.TS;
            obj.X          = R.X;
            obj.comorbid   = R.comorbid;
            obj.pdistr     = R.pdistr;
            obj.totalpdistr = R.totalpdistr;
        end
        function obj = run(obj)
            % normalization (zscore) function handel
            nanz = @(x) (x-nanmean(x))./nanstd(x);
            obj.TS.zage = nanz(str2double(cellstr(obj.TS.age)));  % ZSCORE Age
            % converts numeric label to categorical
            M = containers.Map([1,2],{'DEP1', 'DEP2'});
            obj.TS.Model_Prediction = arrayfun(@(j) M(j),obj.TS.Model_Prediction,'UniformOutput',0);
            % ---------------------
            for parm = ["madrs_total", "bprs_aff", "cains_total", "cains_map", "cains_exp", "tmtA", "tmtB"]
                obj.TS.(strcat(parm,"_delta"))= obj.computedelat(obj.TS, parm);
                % Zscore baseline
                obj.TS.(strcat("z", parm, "_B")) = nanz(obj.TS.(strcat(parm, "_B")));
                % Bootstraping the table
                [obj.mdl.(parm).coeff, obj.mdl.(parm).coeffname, obj.mdl.(parm).idx] = utility.boot_table(obj.TS,...
                    strcat(parm,"_delta ~ 1 + zage + sex + treatment * Model_Prediction +z",parm,"_B"));
                obj.mdl.(parm).eIDX = find(cellfun(@(x) strcmp(x,'treatment_active:Model_Prediction_DEP2'), obj.mdl.(parm).coeffname{1}));
                obj.mdl.(parm).CI  = quantile(obj.mdl.(parm).coeff(:,obj.mdl.(parm).eIDX),[.025,.975]); % 95% bootstrap confidence interval
            end
        end % end run
        function plot(obj)
            iplot = 0;
            clear plot_violn
            for parm = ["madrs_total", "bprs_aff", "cains_total", "cains_map", "cains_exp", "tmtA", "tmtB"]
                iplot = iplot +1;
                % average over bootstrap iterations for ploting
                TP = utility.average_boots_for_plot(obj.mdl,obj.TS,parm,strcat(parm,"_delta"));
                subplot(2,4,iplot)
                cla
                utility.plot_violin(TP,strcat(parm,"_delta"))
                title(strcat("\rm", parm))
            end
        end

    end% end public methods

    methods(Static, Access = private)
        function delta = computedelat(Table, parameter)
            delta = Table.(strcat(parameter,"_F"))- Table.(strcat(parameter,"_B"));
        end
    end % Static and private methods

end