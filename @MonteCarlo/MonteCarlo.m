% -*- coding: 'UTF-8' -*-
classdef MonteCarlo
    % Monte Carlo performs Monte Carlo premutations and probabilistic
    % fitting of the WBM to empirical data

    %   Authors:
    %           Neda Kaboodvand, n.kaboodvand@gmail.com
    %           Behzad Iravani, behzadiravani@gmail.com
    % This function is part of scripts for Macroscopic resting state model predicts
    % theta burst stimulation response: a randomized trial


    properties
        datapath(1,1) cell       % 1x1 cell contains the path to simulated data
        threshold(1,1) double {mustBeLessThan(threshold, 1), mustBeGreaterThan(threshold, 0)} = 0.50;  % a floating number between 0 and 1 that defines the threshold (percentage) randomly selected patients per iteration of Monte Carlo
        numberofperm(1,1) double {mustBeInteger(numberofperm)} = 500;% an integer number that defines the number of permutations
        EachItofSim(1,1) double {mustBeInteger(EachItofSim)}  = 211*68; % a double defines the dimenstion of each simulations timepoints(default, 211) * regions (default, 68)
        randomseed(1,1) double = 15
    end

    properties(Dependent)
        m                      % mapping of the .dat file
        totalS                 % total number of simulations
    end

    methods

        function obj = MonteCarlo(datapath, threshold, numberofperm, EachItoSim, randomseed) % constructor method
            if ~isempty(datapath)
                obj.datapath = datapath;
            else
                error('Please define the path to simulated data.')
            end
            if ~isempty(threshold)
                obj.threshold = threshold;
            end
            if ~isempty(numberofperm)
                obj.numberofperm = numberofperm;
            end
            if ~isempty(EachItoSim)
                obj.EachItofSim = EachItoSim;
            end
            if ~isempty(randomseed)
                obj.randomseed = randomseed;
            end
        end % constructor function
         function run(obj,TCs)
            % Static method of Monte Carlo 
            % Input:
            %       TC  - a cell containing empirical BOLD time series TCs{nsub} 
            runMot(obj, TCs)
        end % end run

        function m = get.m(obj)
            m = memmapfile(obj.datapath{:}, ...
                'Format', 'double' , ...
                'Writable', false);
        end % get method m
        function totalS = get.totalS(obj)
            % get methods for the dependent variable totalS
            totalS = length(obj.m.Data)/obj.EachItofSim;
        end % get method totalS
    end % methods

end % end MonteCarlo