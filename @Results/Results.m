% -*- coding: 'UTF-8' -*-
classdef Results < handle
    % Results produces result figures
    %
    %   Authors:
    %           Neda Kaboodvand, n.kaboodvand@gmail.com
    %           Behzad Iravani, behzadiravani@gmail.com
    % This function is part of scripts for Macroscopic resting state model predicts
    % theta burst stimulation response: a randomized trial
    % Jan 2023, Stanford, Palo Alto, USA


    properties
        MonteCarlo (1,3) cell % a cell containing MonteCarlo permutation data
        pdistr (1,3) cell     % measures of similarity between empirical and simulated data for every threshold
        totalpdistr
        X                     %  simulated timeseries for each sub-cohort
        TS                    %  Table include stratified indexes
        comorbid (1,1) struct  %  chi-squared test for comorbidity and model prediction 
    end % end properties

    methods
        function obj = Results(MonteCarlo)
            obj.MonteCarlo = MonteCarlo;
        end % end constructor


        function obj = plot_parameters_distributions(obj, search)
            unity_normalization = @(x) (x - min(x(:)))/(max(x(:))-min(x(:)));
            p_count = 0;
            close all
            figure('Units', 'Normalized', 'Position', [0 0 1 1])
            for mot = {'30%', '50%', '70%'}
                p_count = p_count +1;
                m = dir(strcat('mot', mot{:}, '\*.mat'));
                for im = 1:length(m)
                    load(strcat('mot', mot{:},'\', m(im).name))

                    [obj.pdistr{p_count}.cor(im), obj.pdistr{p_count}.indx(im)] = max(unity_normalization(corr_staticFC));

                    obj.pdistr{p_count}.FC_corr(im)     = corr_staticFC(obj.pdistr{p_count}.indx(im));
                    obj.pdistr{p_count}.MetaStable(im)  = MSs(obj.pdistr{p_count}.indx(im));
                    obj.pdistr{p_count}.Synchrony(im)   = Sync_s(obj.pdistr{p_count}.indx(im));
                    [i1,i2,i3,i4] = ind2sub([numel(search.A),numel(search.G),numel(search.F),numel(search.M)],obj.pdistr{p_count}.indx(im));
                    obj.pdistr{p_count}.p{im,1} = sprintf('A=%1.3f,G=%1.3f,F=%1.1f,M=%1.2f',...
                        search.A(i1),...
                        search.G(i2),...
                        search.F(i3),...
                        search.M(i4));

                    obj.pdistr{p_count}.Synchrony_emp{im} = mean(Sync_e);
                    obj.pdistr{p_count}.MetaStable_emp{im}= mean(MSe);
                    obj.pdistr{p_count}.FCe{im}           = FCe;

                    obj.pdistr{p_count}.Dis_MS          = abs(MSs-mean(MSe));
                    obj.pdistr{p_count}.Dis_Sync        = abs(Sync_s-mean(Sync_e));
                end
                %% Histogram
                [n,x] = hist(categorical(obj.pdistr{p_count}.p));
                mc = containers.Map(x, 1:numel(x));
                % fit a kernel distribution
                DI = fitdist(cellfun(@(x) mc(x), obj.pdistr{p_count}.p), 'Kernel', 'Kernel', 'normal',  'Width', .25);
                res = 1e3;
                xx  = linspace(1,numel(x),res);
                %%------------------ ploting -----------------------%%
                subplot(1,3,p_count)
                area(xx,...
                    sum(n)*(DI.pdf(xx)))
                [~, mp] = maxk(n,2); % find the peaks
                %%---------------------------------------------------%%
                set(gca, 'XTick', sort(mp), 'XTickLabel', x(sort(mp)),'XTickLabelRotation', 25)
                ax = gca;
                ax.YAxis.Scale = 'line';
                ax.LineWidth     = 1.75;
                ax.TickLength(1) = .02;
                ax.YTick         = ylim();
                ax.YAxis.MinorTick   = 'On';
                ax.YAxis.MinorTickValues = 0:50:600;


                line('XData',[mp(1), mp(1)], 'YData',[0, sum(n)*DI.pdf(mp(1))],'LineWidth',2,'Color', 'w', 'LineStyle', '--')
                line([mp(2), mp(2)], [0, sum(n)*DI.pdf(mp(2))],'LineWidth',2,'Color', 'w', 'LineStyle', '--')

                ylabel('Count')
                title(strcat('\rm Distribution of parameters --', mot{:},'random selection'))
                box off
                axis square
            end
            print result\HistogramParameters.png -dpng

        end % plot paramteres distributions

        function plot_agg_dist(obj)

            %% Histogram
            [n,x] = hist(categorical(obj.totalpdistr.p));
            mc = containers.Map(x, 1:numel(x));
            % fit a kernel distribution
            DI = fitdist(cellfun(@(x) mc(x), obj.totalpdistr.p), 'Kernel', 'Kernel', 'normal',  'Width', .5);
            res = 1e3;
            xx  = linspace(1,numel(x),res);
            %%------------------ ploting -----------------------%%
            area(xx,...
                sum(n)*(DI.pdf(xx)))
            [~, mp] = maxk(n,2); % find the peaks
            %%---------------------------------------------------%%
            set(gca, 'XTick', sort(mp), 'XTickLabel', x(sort(mp)),'XTickLabelRotation', 25)
            ax = gca;
            ax.YAxis.Scale = 'line';
            ax.LineWidth     = 1.75;
            ax.TickLength(1) = .02;
            ax.YTick         = ylim();
            ax.YAxis.MinorTick   = 'On';
            ax.YAxis.MinorTickValues = 0:1e2:2e3;


            line('XData',[mp(1), mp(1)], 'YData',[0, sum(n)*DI.pdf(mp(1))],'LineWidth',2,'Color', 'w', 'LineStyle', '--')
            line([mp(2), mp(2)], [0, sum(n)*DI.pdf(mp(2))],'LineWidth',2,'Color', 'w', 'LineStyle', '--')

            ylabel('Count')
            title(strcat('\rm Distribution of parameters'))
            box off
            axis square
            print result\AggHistogramParameters.png -dpng
            print result\AggHistogramParameters.svg -dsvg
        end % plot_agg_dist

        function obj = plot_synchrony_metastability(obj)
            obj.totalpdistr = concat_(obj.pdistr);
            [n,x] = hist(categorical(obj.totalpdistr.p)); % retrive the histogram for concatenated thresholds
            [np, mp] = maxk(n,2);
            subplot(121)
            clear ylim

            b =  bar(([mean(obj.totalpdistr.MetaStable(categorical(x(mp(1))) == categorical(obj.totalpdistr.p))'), ...
                mean(obj.totalpdistr.MetaStable(categorical(x(mp(2))) == categorical( obj.totalpdistr.p))')]-...
                mean(obj.totalpdistr.MetaStable))./...
                std(obj.totalpdistr.MetaStable));
            b.BarWidth = .5;


            ax = gca;
            ax.LineWidth = 1.75;
            ax.TickLength(1) = .03;
            box off
            axis square
            ylim(1.1*ylim());
            drawnow
            set(gca, 'XTickLabel', {'Dep1', 'Dep2'}, 'YTick', ylim(), 'YTickLabel', round(ylim(),1),...
                'FontName', 'Arial', 'FontSize', 16);
            title('\rm MetaStability')
            ylabel('Normalized change from average')
            xlim([0,3])
            subplot(122)

            b = bar(([mean( obj.totalpdistr.Synchrony(categorical(x(mp(1))) == categorical( obj.totalpdistr.p))'), ...
                mean( obj.totalpdistr.Synchrony(categorical(x(mp(2))) == categorical( obj.totalpdistr.p))')]-...
                mean( obj.totalpdistr.Synchrony))./...
                std( obj.totalpdistr.Synchrony));% % DIS = fitdist( obj.totalpdistr.Synchrony', 'Kernel', 'Kernel', 'normal',  'Width', 0.003);

            b.BarWidth = .5;

            ax = gca;

            ax.LineWidth = 1.75;
            ax.TickLength(1) = .03;
            box off
            axis square
            ylim(1.1*ylim());
            set(gca, 'XTickLabel', '' , 'YTick', ylim(), 'YTickLabel', round(ylim(),2),...
                'FontName', 'Arial', 'FontSize', 16);
            title('\rm Synchrony')
            xlim([0,3])
            % ylabel('Percentage change from average')
            print result\barMeasures.svg -dsvg
            print result\barMeasures.png -dpng

            % spider ploting parameters
            close all
            SpiderPlot_parameters([cellfun(@str2num,regexp(x{mp(1)},'(-?)\d.\d+','match'))
                cellfun(@str2num,regexp(x{mp(2)},'(-?)\d.\d+','match'))]);


            function mxc = concat_(mx)
                % contcatenate the measures
                mx0 = cat(1, mx{:});
                mx = mx0;
                clear mx0

                for f = fieldnames(mx)'
                    [s0,~] = size(mx(1).(f{:}));
                    if s0==1, mxc.(f{:}) = cat(2,mx.(f{:}));
                    else,mxc.(f{:}) = cat(1,mx.(f{:}));
                    end
                end
            end % concat


            function SpiderPlot_parameters(P)

                addpath(genpath('C:\MatlabToolboxes\spider'))
                P(:,3) = abs(P(:,3)); % The sign has incoporated in the formula.
                % DEP1{'A=0.010, G = 0.020, f = -0.200, M = 0.110'}
                % DEP2{'A=0.054, G = 0.026, f = -0.300, M = 0.190'}

                % P  = [0.01, .02, .2, 0.11;0.054, .026, .3, 0.19];
                spider_plot(P,...
                    'AxesLabels', {'A', 'G', 'F', 'M'},...
                    'AxesInterval', 2,...
                    'AxesLimits',[-0.2, 0.01, 0.2, 0.0; 0.2 0.04 1 0.5],...
                    'AxesPrecision',2,...
                    'FillOption',{'on','on'}, 'Color', [0 77 129;237 137 93]./255);

                print -dsvg result\optimizedPar.svg
                print -dpng result\optimizedPar.png

            end % spider
        end % plot
        function obj = Stratifying(obj, Empirical, T)

            [n,x] = hist(categorical(obj.totalpdistr.p)); % retrive the histogram for concatenated thresholds
            [~, mp] = maxk(n,2);


            nS = unique([obj.totalpdistr.indx(categorical(x(mp(1))) == categorical(obj.totalpdistr.p))]);

            obj.X.dep1.TC = reshape(...
                obj.MonteCarlo{1}.m.Data((1+(nS-1)*obj.MonteCarlo{1}.EachItofSim):((nS)*obj.MonteCarlo{1}.EachItofSim)),...
                211,68);

            nS = unique([obj.totalpdistr.indx(categorical(x(mp(2))) == categorical(obj.totalpdistr.p))]);

            obj.X.dep2.TC = reshape(...
                obj.MonteCarlo{1}.m.Data((1+(nS-1)*obj.MonteCarlo{1}.EachItofSim):((nS)*obj.MonteCarlo{1}.EachItofSim)),...
                211,68);

            [obj.X.dep1.M, obj.X.dep1.S] = simple_comp_MS(obj.X.dep1.TC);
            [obj.X.dep2.M, obj.X.dep2.S] = simple_comp_MS(obj.X.dep2.TC);

            obj.X.dep1.FC = simple_comp_FC(obj.X.dep1.TC);
            obj.X.dep2.FC = simple_comp_FC(obj.X.dep2.TC);
            FCe = mean(cat(3, obj.totalpdistr.FCe{:}),3);
            obj.X.dep1.sim = corr(squareform(FCe)',...
                squareform(obj.X.dep1.FC)');
            obj.X.dep2.sim =corr(squareform(FCe)',...
                squareform(obj.X.dep2.FC)');
            clear E
            for ie = 1:size(Empirical.TC,2)

                [E(ie).M, E(ie).S ] = simple_comp_MS(Empirical.TC{ie});

                DIST(:,1,ie) = [1-(E(ie).M-obj.X.dep1.M)/max(obj.totalpdistr.Dis_MS)
                    1-(E(ie).S-obj.X.dep1.S)/max(obj.totalpdistr.Dis_Sync)
                    corr(squareform(simple_comp_FC(Empirical.TC{ie}))',...
                    squareform(obj.X.dep1.FC)')];

                DIST(:,2,ie) = [1-(E(ie).M-obj.X.dep2.M)/max(obj.totalpdistr.Dis_MS)
                    1-(E(ie).S-obj.X.dep2.S)/max(obj.totalpdistr.Dis_Sync)
                    corr(squareform(simple_comp_FC(Empirical.TC{ie}))',...
                    squareform(obj.X.dep2.FC)')];
            end

            [~,cid]= max(squeeze(DIST(3,:,:)),[],1); % used only FC dimension 3

            % find DEP subjects
            obj.TS = T(~isnan(T.fMRIses1) & categorical(cellstr(T.group)) == 'depression',:);

            obj.TS.Model_Prediction  = cid';
            obj.TS.MetaStability     = [E.M]';
            obj.TS.Synchrony         = [E.S]';
            writetable(obj.TS,'result\TabularData_Sep25_2022.csv')

        end % Stratified

        function anova_movement_modelpred()


        end % anova_movement_modelpred

        function obj = chi2_comorbid(obj, comorbid)
            obj.TS = obj.add_comorbid(comorbid, obj.TS);
            %% Anxiety
            [obj.comorbid.table.Anxiety, obj.comorbid.chi2.Anxiety, obj.comorbid.p.Anxiety, ~] = crosstab(obj.TS.ComorbidityAnxietyDisorder, obj.TS.Model_Prediction);
            %% Austism
            [obj.comorbid.table.Austism, obj.comorbid.chi2.Austism, obj.comorbid.p.Austism, ~] = crosstab(obj.TS.ComorbidityAutismSpectrumDisorder, obj.TS.Model_Prediction);
            %% ADHD
            [obj.comorbid.table.ADHD, obj.comorbid.chi2.ADHD, obj.comorbid.p.ADHD, ~] = crosstab(obj.TS.ComorbidityADHD_ADD, obj.TS.Model_Prediction);
            %% PersonalityDisorder
            [obj.comorbid.table.PersonalityDisorder, obj.comorbid.chi2.PersonalityDisorder, obj.comorbid.p.PersonalityDisorder, ~] = crosstab(obj.TS.ComorbidityPersonalityDisorder, obj.TS.Model_Prediction);
            %% Bipolar
            [obj.comorbid.table.Bipolar, obj.comorbid.chi2.Bipolar, obj.comorbid.p.Bipolar, ~] = crosstab(obj.TS.BipolarDisorder, obj.TS.Model_Prediction);     
        end% chi2_comorbid test

    end % methods


end