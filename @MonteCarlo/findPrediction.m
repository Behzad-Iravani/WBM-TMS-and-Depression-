clc
clear all
close all
% search space

search = struct('A', single([-.2:.01:.2]), 'G',single([0.005:0.005:0.05]),...
    'F', single([-1:0.1:0]), 'M', single([0.1:0.03:.5]));


% search = struct('A', -0.044:0.002:0.044, 'G', 0.002:0.002:0.026, 'F', -1:0.1:0, 'M', 0.25:0.01:0.42);
% normlization function
unity_normalization = @(x) (x - min(x(:)))/(max(x(:))-min(x(:)));

% % % old_i(1) = sub2ind([numel(search.A),numel(search.G),numel(search.F),numel(search.M)],...
% % %     find(search.A == -0.04), find(search.G == .014), find(search.F == -.4), find(search.M == .30))
% % % old_i(2) = sub2ind([numel(search.A),numel(search.G),numel(search.F),numel(search.M)],...
% % %     find(search.A == -0.02), find(abs(search.G - 0.020000000000000)<=eps), find(search.F == -.4), find(search.M == .32))
% % % % corr_staticFC(old_i)
p_count = 0;
close all
figure('Units', 'Normalized', 'Position', [0 0 1 1])
for mot = {'30%', '50%', '70%'}
    p_count = p_count +1;
    m = dir(strcat('mot', mot{:}, '\*.mat'));
    for im = 1:length(m)
        load(strcat('mot', mot{:},'\', m(im).name))

        [mx{p_count}.cor(im), mx{p_count}.indx(im)] = max(unity_normalization(corr_staticFC));

        mx{p_count}.FC_corr(im)     = corr_staticFC(mx{p_count}.indx(im));
        mx{p_count}.MetaStable(im)  = MSs(mx{p_count}.indx(im));
        mx{p_count}.Synchrony(im)   = Sync_s(mx{p_count}.indx(im));
        [i1,i2,i3,i4] = ind2sub([numel(search.A),numel(search.G),numel(search.F),numel(search.M)],mx{p_count}.indx(im));
        mx{p_count}.p{im,1} = sprintf('A=%1.3f,G=%1.3f,F=%1.1f,M=%1.2f',...
            search.A(i1),...
            search.G(i2),...
            search.F(i3),...
            search.M(i4));

        %     reshape(corr_staticFC, numel(search.A),numel(search.G),numel(search.F),numel(search.M))
        mx{p_count}.Synchrony_emp{im} = mean(Sync_e);
        mx{p_count}.MetaStable_emp{im}= mean(MSe);
        mx{p_count}.FCe{im}           = FCe;
    end
    %% Histogram
    [n,x] = hist(categorical(mx{p_count}.p));
    mc = containers.Map(x, 1:numel(x));
    % histogram(categorical(mx{p_count}.p))
    DI = fitdist(cellfun(@(x) mc(x), mx{p_count}.p), 'Kernel', 'Kernel', 'normal',  'Width', .25);

    res = 1e3;
    xx  = linspace(1,numel(x),res);
    % pp  = interp1(1:numel(x),DI.pdf(1:numel(x)), xx,'pchip');
    % pp  = (pp/(sum(pp)))*sum(n);
    sum(n)
    sum(DI.pdf(xx))
    subplot(1,3,p_count)
    area(xx,...
        sum(n)*(DI.pdf(xx)))


    [~, mp] = maxk(n,2);
    % %
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
print HistogramParameters.svg -dsvg

nanmean(mx.FC_corr(categorical(x(mp(1))) == categorical(mx.p)))
nanmean(mx.FC_corr(categorical(x(mp(2))) == categorical(mx.p)))

quantile([mx.Synchrony_emp{categorical(x(mp(1))) == categorical(mx.p)}],[.5])
quantile([mx.Synchrony_emp{categorical(x(mp(2))) == categorical(mx.p)}],[.5])

quantile([mx.MetaStable_emp{categorical(x(mp(1))) == categorical(mx.p)}],[.5])
quantile([mx.MetaStable_emp{categorical(x(mp(2))) == categorical(mx.p)}],[.5])




%%

% subplot(131)
% boxplot([mx.FC_corr(categorical(x(mp(1))) == categorical(mx.p))'; ...
%     mx.FC_corr(categorical(x(mp(2))) == categorical(mx.p))'], [ones(sum(categorical(x(mp(1))) == categorical(mx.p)),1)
%     2*ones(sum(categorical(x(mp(2))) == categorical(mx.p)),1)],...
%     'BoxStyle','filled','MedianStyle','target','Colors','k','Symbol',' ')
% % DIFC = fitdist(mx.FC_corr', 'Kernel', 'Kernel', 'normal',  'Width', 0.008);
% % area(linspace(.15,.3,1e2),...
% %     DIFC.pdf(linspace(.15,.3,1e2))...
% %     )
% ax = gca;
% % ax.YAxis.Scale = 'log';
% ax.LineWidth = 1.75;
% ax.TickLength(1) = .03;
% box off
% axis square
% set(gca, 'XTickLabel', {'Dep1', 'Dep2'}, 'YTick', ylim(), 'YTickLabel', round(ylim(),2))
subplot(121)
clear ylim
% boxplot([[mx.MetaStable_emp{categorical(x(mp(1))) == categorical(mx.p)}]'; ...
%     [mx.MetaStable_emp{categorical(x(mp(2))) == categorical(mx.p)}]'], [ones(sum(categorical(x(mp(1))) == categorical(mx.p)),1)
%     2*ones(sum(categorical(x(mp(2))) == categorical(mx.p)),1)],...
%     'BoxStyle','filled','MedianStyle','target','Colors','k','Symbol',' ')

b =  bar(1e2*([mean(mx.MetaStable(categorical(x(mp(1))) == categorical(mx.p))'), ...
    mean(mx.MetaStable(categorical(x(mp(2))) == categorical(mx.p))')]-...
    mean(mx.MetaStable))./...
    mean(mx.MetaStable))
b.BarWidth = .5;

% % DIM = fitdist(mx.MetaStable', 'Kernel', 'Kernel', 'normal',  'Width', 0.001);
% %
% % area(linspace(0.1,.15,1e2),...
% %     DIM.pdf(linspace(0.1,.15,1e2))...
% %     )
ax = gca;
ax.LineWidth = 1.75;
ax.TickLength(1) = .03;
box off
axis square
ylim(1.1*ylim());
drawnow
set(gca, 'XTickLabel', {'Dep1', 'Dep2'}, 'YTick', ylim(), 'YTickLabel', round(ylim(),1))
title('\rm MetaStability')
ylabel('Percentage change from average')

subplot(122)
% boxplot([[mx.Synchrony_emp{categorical(x(mp(1))) == categorical(mx.p)}]'; ...
%     [mx.Synchrony_emp{categorical(x(mp(2))) == categorical(mx.p)}]'], [ones(sum(categorical(x(mp(1))) == categorical(mx.p)),1)
%     2*ones(sum(categorical(x(mp(2))) == categorical(mx.p)),1)],...
%     'BoxStyle','filled','MedianStyle','target','Colors','k','Symbol',' ')

b = bar(1e2*([mean(mx.Synchrony(categorical(x(mp(1))) == categorical(mx.p))'), ...
    mean(mx.Synchrony(categorical(x(mp(2))) == categorical(mx.p))')]-...
    mean(mx.Synchrony))./...
    mean(mx.Synchrony));% % DIS = fitdist(mx.Synchrony', 'Kernel', 'Kernel', 'normal',  'Width', 0.003);

b.BarWidth = .5;
% % area(linspace(0.2,.3,1e2),...
% %     DIS.pdf(linspace(0.2,.3,1e2))...
% %     )
ax = gca;
% ax.YAxis.Scale = 'log';
ax.LineWidth = 1.75;
ax.TickLength(1) = .03;
box off
axis square
ylim(1.1*ylim());
set(gca, 'XTickLabel', '' , 'YTick', ylim(), 'YTickLabel', round(ylim(),2))
title('\rm Synchrony')
% ylabel('Percentage change from average')


print HistogramMeasures.svg -dsvg
%%

addpath(genpath('C:\MatlabToolboxes\spider'))
P = [cellfun(@str2num,regexp(x{mp(1)},'(-?)\d.\d+','match'))
    cellfun(@str2num,regexp(x{mp(2)},'(-?)\d.\d+','match'))];
optimized_paramertes_plot(P)

% ind[find(abs(search.A - P(1,1))<1e-3), find(search.G == P(1,2)), find(search.F == -P(1,3)), find(search.M == P(1,4))]

% histogram(categorical(p),Normalization="probability")
% % % % cx =0;
% % % % clear dd
% % % % for xx = x
% % % %     cx =cx+1;
% % % %     ind  = cellfun(@(s) strcmp(s,xx), p);
% % % %     if sum(ind)~=1
% % % %         dd(cx,:) = mean(d(ind,:));
% % % %     else
% % % %         dd(cx,:) = d(ind,:);
% % % %     end
% % % % end
% % % %
% % % % squareform(pdist(dd))
% % % % bar(n)
% % % % set(gca,'XTick',1:15, 'XTickLabel', x)
% % % % % fitted parameters
% % % % ;
%%
% % % % % model = WBM;
% % % % % model.intiate_model
% % % % %
% % % % % load Data\T_plusBehavior.mat
% % % % % % find DEP subjects
% % % % % T = T(~isnan(T.fMRIses1) & categorical(cellstr(T.group)) == 'depression',:);
% % % % %
% % % % % subjs = dir('TCs_subjDK\ses-baseline\subj*');
% % % % % c = 0;
% % % % % for sub = T.ID'
% % % % %     c = c + 1;
% % % % %     fprintf('load TCs for %s\n', sub{:})
% % % % %     %     load(strcat("data\",sub{:},filesep, "Atlas_wmparcTCs_ses1.mat"))
% % % % %     ind = find(cellfun(@(x) contains(x,sub{:}), {subjs.name}));
% % % % %     load(fullfile(subjs(ind).folder, subjs(ind).name, 'TCs.mat'))
% % % % %     load(fullfile(subjs(ind).folder, subjs(ind).name, 'ROIs.mat'))
% % % % %
% % % % %     % extract cortical time courses
% % % % %     TC{c}   = utility.find_cortical_nodes(TCs, ROIs);
% % % % %     FQ(c,:) = utility.find_freq(TC{c});
% % % % % end
% % % % % disp('Average the empirical data')
% % % % % Empirical.treatment = T.treatment;
% % % % % Empirical.TC    =  TC;
% % % % % Empirical.Conn  = arrayfun(@(j) corrcoef(Empirical.TC{j}),1:length(Empirical.TC),'UniformOutput',false);
% % % % %
% % % % % % Empirical.Conn  = simple_comp_FC(mean(cat(3, Empirical.Conn{:}),3));
% % % % %
% % % % % % frequecny of BOLD activity
% % % % % Empirical.freq  = median(FQ);
% % % % % % fsample
% % % % % Empirical.fsamp = 1/2; % 1 over TR
% % % % %
% % % % % fnyq   = Empirical.fsamp/2;
% % % % % [b, a] = butter(model.filter.order, [model.filter.lp, model.filter.hp]./fnyq);
% % % % %
% % % % %
% % % % % % add empirical data
% % % % % model.Empirical = Empirical;
% % % % %
% % % % % model.parameters.A = P(1,1);
% % % % % model.parameters.G = P(1,2);
% % % % % model.parameters.F = -P(1,3); % the negative sign is needed.
% % % % % model.parameters.M = P(1,4);
% % % % % clear X Y w
% % % % % [X.dep1,Y.dep1,w.dep1] = model.de_Simulate_tired2_v22_nodalInput();
% % % % % X.dep1    = filtfilt(b,a,X.dep1);
% % % % %
% % % % %
% % % % % model.parameters.A = P(2,1);
% % % % % model.parameters.G = P(2,2);
% % % % % model.parameters.F = -P(2,3); % the negative sign is needed.
% % % % % model.parameters.M = P(2,4);
% % % % %
% % % % % [X.dep2,Y.dep2,w.dep2] = model.de_Simulate_tired2_v22_nodalInput();
% % % % % X.dep2    = filtfilt(b,a,X.dep2);
% % % % %
% % % % % subplot(221)
% % % % % plot(conv(max(xai(:,:,abs(search.F - -P(1,3))<1e-3, search.M == P(1,4),3),[],2),ones(1,5).*.2,'same'))
% % % % % % axis off
% % % % % pbaspect([1,.25,1])
% % % % % subplot(222)
% % % % % plot(conv(max(xai(:,:,search.F == -P(2,3), search.M == P(2,4),3),[],2),ones(1,5).*.2,'same'))
% % % % % % axis off
% % % % % pbaspect([1,.25,1])
% % % % % subplot(223)
% % % % % contourf(search.G,search.A,...
% % % % %     conv2(squeeze(xa(:,:, search.F == -P(1,3), search.M == P(1,4))),...
% % % % %     [1 1;1 1]/4, 'same'),4, '--')
% % % % % hold on
% % % % % scatter(P(1,2),P(1,1), 120,'x', 'LineWidth', 2, 'MarkerEdgeColor', 'r')
% % % % % subplot(224)
% % % % % contourf(search.G,search.A,...
% % % % %     conv2(squeeze(xa(:,:, search.F == -P(2,3), search.M == P(2,4))),...
% % % % %     [1 1 ;1 1]/4, 'same'),4, '--')
% % % % % hold on
% % % % % scatter(P(2,2),P(2,1), 120,'x', 'LineWidth', 2, 'MarkerEdgeColor', 'r')
% % % % %
% % % % %
% % % % % % % [mE, sE]= cellfun(@(x) simple_comp_MS(x), Empirical.TC);
% % % % % % % mean(mE)
% % % % % % % mean(sE)

mf = memmapfile('Data\simulated.dat', ...
    'Format', 'double' , ...
    'Writable', false);
EachItofSim = 211*68;

nS = unique([mx.indx(categorical(x(mp(1))) == categorical(mx.p))]);

X.dep1 = reshape(...
    mf.Data((1+(nS-1)*EachItofSim):((nS)*EachItofSim)),...
    211,68);

nS =  unique([mx.indx(categorical(x(mp(2))) == categorical(mx.p))]);

X.dep2 = reshape(...
    mf.Data((1+(nS-1)*EachItofSim):((nS)*EachItofSim)),...
    211,68);

[d1.M, d1.S] = simple_comp_MS(X.dep1);
[d2.M, d2.S] = simple_comp_MS(X.dep2);
d1.FC = simple_comp_FC(X.dep1);
d2.FC = simple_comp_FC(X.dep2);
corr(squareform(FCe)',...
    squareform(d1.FC)')

corr(squareform(FCe)',...
    squareform(d2.FC)')



load Data\Empirical.mat
for ie = 1:size(Empirical.TC,2)
    corr(squareform(simple_comp_FC(Empirical.TC{ie}))',...
        squareform(d1.FC)')
    clear E
    [E.M, E.S ] = simple_comp_MS(Empirical.TC{ie});

    DIST(:,1,ie) = [1-(E.M-d1.M)/max(Dis_MS)
        1-(E.S-d1.S)/max(Dis_Sync)
        corr(squareform(simple_comp_FC(Empirical.TC{ie}))',...
        squareform(d1.FC)')];

    DIST(:,2,ie) = [1-(E.M-d2.M)/max(Dis_MS)
        1-(E.S-d2.S)/max(Dis_Sync)
        corr(squareform(simple_comp_FC(Empirical.TC{ie}))',...
        squareform(d2.FC)')];
end

[~,cid]= max(squeeze(DIST(3,:,:)),[],1);
load Data\T_plusBehavior.mat
% find DEP subjects
T = T(~isnan(T.fMRIses1) & categorical(cellstr(T.group)) == 'depression',:);

T.Model_Prediction  = cid';
writetable(T,'TabularData_Sep25_2022.csv')
% %
% % subplot(2,4,[2,3])
% % imagesc(Empirical.Conn)
% % axis square
% % axis off
% % subplot(2,4,[5,6])
% % imagesc(simple_comp_FC(X.dep1))
% % axis square
% % axis off
% % subplot(2,4,[7,8])
% % imagesc(simple_comp_FC(X.dep2))
% % axis square
% % axis off
% % corr(squareform(Empirical.Conn)',squareform(simple_comp_FC(X.dep1))')
% % corr(squareform(Empirical.Conn)',squareform(simple_comp_FC(X.dep2))')
% % corr(squareform(simple_comp_FC(X.dep1))',squareform(simple_comp_FC(X.dep2))')
%%
% col = lines(2);
% subplot(121)
% ph = polarhistogram(mean(angle(hilbert(X.dep1'))));
% ph.FaceColor = col(1,:);
% ph.FaceAlpha = .85;
% title('\rmDEP1')
% subplot(122)
% ph = polarhistogram(mean(angle(hilbert(X.dep2'))));
% ph.FaceColor = col(2,:);
% ph.FaceAlpha = .85;
% title('\rmDEP2')
%
% [h, p] = kstest2(mean(angle(hilbert(X.dep1'))),mean(angle(hilbert(X.dep2'))))
