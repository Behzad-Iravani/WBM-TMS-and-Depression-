% -*- coding: 'UTF-8' -*-
function [MS, Sync] = simple_comp_MS(TC)
% simple_comp_MS is a private method of Results that computes synchrony and 
% metastability of time course TC using Hilbert transform
%   Authors:
%           Neda Kaboodvand, n.kaboodvand@gmail.com
%           Behzad Iravani, behzadiravani@gmail.com
% This function is part of scripts for Macroscopic resting state model predicts
% theta burst stimulation response: a randomized trial

IP = angle(hilbert(TC))';%Instantaneous Phases of NxT fMRI data: node(N) by time(T).
%If Xr is a matrix, then HILBERT operates along the columns of Xr.
MS = std(abs(sum(exp(1i.*IP),1)./size(IP,1)));
Sync = mean(abs(sum(exp(1i.*IP),1)./size(IP,1)));
end