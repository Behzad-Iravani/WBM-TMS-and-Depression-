% -*- coding: 'UTF-8' -*-
function runMot(obj,TCs)
% runMot is a private method of MonteCarlo that performs ther permutations
%
%   Authors:
%           Neda Kaboodvand, n.kaboodvand@gmail.com 
%           Behzad Iravani, behzadiravani@gmail.com 
% This function is part of scripts for Macroscopic resting state model predicts
% theta burst stimulation response: a randomized trial

%% Monte Carlo
rng(obj.randomseed)
% Mot.randindex   = randi(numel(Empirical.TC),Mot.number_iter,floor(numel(Empirical.TC)/2));
Mot.randindex   = rand(obj.numberofperm,numel(TCs))>(1-obj.threshold);
%% ----------------------------------- MonteCarlo -----------------------------------------------------------
D = parallel.pool.DataQueue;
afterEach(D, @(it) fprintf('%d iteration out of %d. Elapsed time %1.3f \n', it(1), it(2), it(3)))

% check for existing runs
prev_runs = dir(strcat(sprintf('mot%d%%', fix(1e2*obj.threshold)),'\measure*.mat'));
nPrev_run = numel(prev_runs);
fprintf('%d runs has been found.\n',nPrev_run)
for mot_ = 1:obj.numberofperm
    t = tic();
    mo = matfile([sprintf('mot%d%%', fix(1e2*obj.threshold)), filesep, sprintf('measures%d', mot_+nPrev_run), '.mat'],'Writable',true);
    FC_e = comp_FC(TCs, Mot.randindex(mot_,:),'E'); % compute empirical FC
    [MS_e, Sync_e] = comp_MS(TCs, Mot.randindex(mot_,:),'E'); % compute empirical MS


    corr_staticFC = nan(obj.totalS ,1);
    MS_s          = nan(obj.totalS ,1);
    Sync_s        = nan(obj.totalS ,1);

    parfor nS = 1:obj.totalS  % total number of simulation
        FC_s = comp_FC(reshape(...
            obj.m.Data((1+(nS-1)*obj.EachItofSim):((nS)*obj.EachItofSim)),...
            211,68),...
            [], 'S'); % compute simulate FC
        if ~all(isnan(FC_s))
            corr_staticFC(nS) =  corr(squareform(FC_e)',squareform(FC_s)');
            
            [MS_s(nS), Sync_s(nS)] = comp_MS(reshape(...
                obj.m.Data((1+(nS-1)*obj.EachItofSim):((nS)*obj.EachItofSim)),...
                211,68),...
                [], 'S'); % compute simulated MS
        end
    end
    mo.corr_staticFC = corr_staticFC;
    % Empirical
    mo.FCe           = FC_e;
    mo.MSe           = MS_e;
    mo.Sync_e        = Sync_e;
    % Simulation
    mo.MSs           = MS_s;
    mo.Sync_s        = Sync_s;
%     mo.FCs           = FC_s;

    mo.Dis_MS          = abs(MS_s-mean(MS_e));
    mo.Dis_Sync        = abs(Sync_s-mean(Sync_e));

    send(D,[mot_, obj.numberofperm,   toc(t)])

end
end

function FCo = comp_FC(TC, index, mod)
Z_fc = @(CC) .5*[log(.00000001+1+CC)-log(.00000001+1-CC)]; % fisher z transformation
if strcmp(mod, 'E')
    TC_tmp = cat(3,TC{index});
    for sub = 1:size(TC_tmp,3)
        FC = corrcoef(squeeze(TC_tmp(:,:,sub)));%210x68
        FC= (~eye(size(FC,1))).*FC; % zero the diagonal element
        FC = Z_fc(FC);
        FC(FC<0)=0; % remove negative FCs
        FC_All(:,:,sub) = FC;
    end
    FCo = mean(FC_All,3);
elseif strcmp(mod, 'S')
    FC = corrcoef(TC);%210x68
    FC= (~eye(size(FC,1))).*FC; % zero the diagonal element
    FCo = Z_fc(FC);
    FCo(FCo<0)=0; % remove negative FCs
end
end


function [MS, Sync] = comp_MS(TC, index, mod)
if strcmp(mod, 'E')
    TC_tmp = cat(3,TC{index});
    for sub = 1:size(TC_tmp,3)
        IP = angle(hilbert(squeeze(TC_tmp(:,:,sub))))';%Instantaneous Phases of NxT fMRI data: node(N) by time(T).
        %If Xr is a matrix, then HILBERT operates along the columns of Xr.
        MS(:,sub) = std(abs(sum(exp(1i.*IP),1)./size(IP,1)));
        Sync(:,sub) = mean(abs(sum(exp(1i.*IP),1)./size(IP,1)));
    end
else
    IP = angle(hilbert(TC))';%Instantaneous Phases of NxT fMRI data: node(N) by time(T).
    %If Xr is a matrix, then HILBERT operates along the columns of Xr.
    MS = std(abs(sum(exp(1i.*IP),1)./size(IP,1)));
    Sync = mean(abs(sum(exp(1i.*IP),1)./size(IP,1)));
end

end
