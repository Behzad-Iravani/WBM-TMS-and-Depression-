% -*- coding: 'UTF-8' -*-
function FCo = simple_comp_FC(TC)
% simple_comp_FC is a private method of Results that computes simple FC of
% time course TC using prearson correlation.
%   Authors:
%           Neda Kaboodvand, n.kaboodvand@gmail.com
%           Behzad Iravani, behzadiravani@gmail.com
% This function is part of scripts for Macroscopic resting state model predicts
% theta burst stimulation response: a randomized trial
Z_fc = @(CC) .5*[log(.00000001+1+CC)-log(.00000001+1-CC)]; % fisher z transformation
[n1,n2] = size(TC);
if n1 ~= n2
    FC = corrcoef(TC);%210x68
else
    FC = TC;
end
FC= (~eye(size(FC,1))).*FC; % zero the diagonal element
FCo = Z_fc(FC);
FCo(FCo<0)=0;
end