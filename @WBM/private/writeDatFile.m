% -*- coding: 'UTF-8' -*-
function writeDatFile(path,Simulated)
% writeDatFile is private method of WBM that writes the simulted data to a
% .dat file. 
%   Inputs:
%           path        -the path to .dat file wherein the simulation is stored
%           Simulated   -the structure that contains simulation.
%   Authors:
%           Neda Kaboodvand, n.kaboodvand@gmail.com 
%           Behzad Iravani, behzadiravani@gmail.com 
% This function is part of scripts for Macroscopic resting state model predicts
% theta burst stimulation response: a randomized trial

%%--------------------------------------------------------------%%
dat = cat(3,Simulated.X{:});
fid = fopen(path,'w');
fwrite(fid, dat, 'double')
fclose(fid)
end