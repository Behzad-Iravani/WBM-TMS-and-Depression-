% -*- coding: 'UTF-8' -*-
classdef WBM < handle
    % WBM performs dynamic system modeling using adaptive frequency hopf model
    % introduced by Kaboodvand et al. Network science 2019,and heun integration method.
    % 
    % 
    % WBM properties:
    %  parameters - the Adaptive frequency hopf model has 4 paramaters
    %               including, A, global bifurcation, G, global coupling,
    %               F, feedback and M, modulation.  
    %  search     - determines the serach range of parameter space 
    %  random     - defines the random seed for reproducibility and the standar deviation of the 
    %               Wiener process (beta, default value is 0.02)
    %  dt         - time step, default is 0.1 (s) 
    %  SC         - structural connectivity
    %  time       - defines the the duration and steps of simulation 
    %  Emiprical  - contains the empirical data 
    %  Simulated  - stores the simualtion data
    %  filter     - temporal filter specification
    % .
    % 
    %   Authors:
    %           Neda Kaboodvand, n.kaboodvand@gmail.com 
    %           Behzad Iravani, behzadiravani@gmail.com 
    %
    % Reference: Neda Kaboodvand, Martijn P. van den Heuvel, Peter Fransson; 
    % Adaptive frequency-based modeling of whole-brain oscillations:
    % Predicting regional vulnerability and hazardousness rates. Network Neuroscience 2019;
    % 3 (4): 1094–1120. doi: https://doi.org/10.1162/netn_a_00104

    properties
        % model's parameter
        parameters = struct('A', [], 'G', [], 'F', [], 'M',[])
        % search range
        search = struct('A', -0.07:0.002:0.07, 'G', 0.002:0.002:0.03, 'F', -0.6:0.1:-0.2, 'M', 0.1:0.01:0.2)
        % modeling random processes
        random     = struct('seed', 15, 'RR', [], 'beta', 0.02)
        % time resolution
        dt = 0.1;
        % structual connectivity
        SC = [];
        % number of regions to simulate
        num_nodes = 68;
        % Time parameters
        time     = struct('dur', 8*60, 't_0', 0, 't_steps', []);
        % empirical data
        Empirical     = struct('TC', [], 'Conn', [], 'freq',[] , 'fsamp', [], 'treatment', []);
        % simulated data
        Simulated     = struct('X', [], 'Conn', [])
        % bandpass filtering
        filter        = struct('lp',.008, 'hp', .1, 'order', 2)

    end

    methods
        function obj = intiate_model(obj)
            % This function intitates model
            % set the time parameters
            % duration of simulation
            obj.time.t_steps = round((obj.time.dur-obj.time.t_0)/obj.dt+1);
            % loading structual connectivity
            load @WBM\+Be\SuperCon68.mat SuperCon3
            obj.SC = SuperCon3;
            % The Wiener process that models the ranodom effects in the brain
            rng(obj.random.seed)
            obj.random.RR = randn(obj.num_nodes,obj.time.t_steps);


        end
        function obj = simulate(obj)
            % simulate is method of class WBM and carries out the
            % simulation
            % setting up data queue
            D = parallel.pool.DataQueue;
            it_total = length(obj.search.A)*length(obj.search.G)*...
                length(obj.search.F)*length(obj.search.M);
            fprintf('total iteration is %d\n\n', it_total)
            afterEach(D, @update_n)
            % --- function handles
            handle_sim    = @(obj, A, G, F, M) sim(obj, A, G, F, M);
            handle_filter = @(obj,X) bps_filter(obj,X);
            % ---
            ci = 0;
            loop_c = 0;
            Xout = cell(1,1);
            Yout = cell(1,1);
            wout = cell(1,1);
            parfor A = 1:length(obj.search.A)%for A = 1:length(obj.search.A)%
                tmpxG = cell(1,length(obj.search.G));
                tmpyG = cell(1,length(obj.search.G));
                tmpwG = cell(1,length(obj.search.G));
                for G = 1:length(obj.search.G)
                    tmpxF = cell(1,length(obj.search.F));
                    tmpyF = cell(1,length(obj.search.F));
                    tmpwF = cell(1,length(obj.search.F));
                    for F = 1:length(obj.search.F)
                        tmpxM = cell(1,length(obj.search.M));
                        tmpyM = cell(1,length(obj.search.M));
                        tmpwM = cell(1,length(obj.search.M));
                        for M = 1:length(obj.search.M)
                            ci = ci +1;
                            % ---
                            [X,Y,w] = ...
                                handle_sim(obj, A, G, F, M);
                            % filtering
                            if ~any(isnan(X))
                                X = handle_filter(obj, X);
                            end
                            if ~any(isnan(Y))
                                Y = handle_filter(obj, Y);
                            end
                            tmpxM{M} = X;
                            tmpyM{M} = Y;
                            tmpwM{M} = w;
                            send(D,A);
                        end
                        tmpxF{F} = tmpxM;
                        tmpyF{F} = tmpyM;
                        tmpwF{F} = tmpwM;
                    end
                    tmpxG{G} = tmpxF;
                    tmpyG{G} = tmpyF;
                    tmpwG{G} = tmpwF;
                end
                Xout{A} = tmpxG;
                Yout{A} = tmpyG;
                wout{A} = tmpwG;

            end


            %% unnesting cells
            while any(cellfun(@iscell,Xout))
                Xout = [Xout{cellfun(@iscell,Xout)} Xout(~cellfun(@iscell,Xout))];
            end
            while any(cellfun(@iscell,Yout))
                Yout = [Yout{cellfun(@iscell,Yout)} Yout(~cellfun(@iscell,Yout))];
            end
            while any(cellfun(@iscell,wout))
                wout = [wout{cellfun(@iscell,wout)} wout(~cellfun(@iscell,wout))];
            end
            obj.Simulated.X = Xout;
            obj.Simulated.Y = Yout;
            obj.Simulated.w = wout;

            % write simulation to file
            writeDatFile('Data\simulated.dat', obj.Simulated)

           
            function [X,Y,w] = sim(obj,A, G, F, M)
                % --- update parameters
                obj.parameters.A = obj.search.A(A);
                obj.parameters.G = obj.search.G(G);
                obj.parameters.F = obj.search.F(F);
                obj.parameters.M = obj.search.M(M);
                [X,Y,w] = ...
                    de_Simulate_v22_nodalInput(obj);
            end
            function update_n(~)
                loop_c = loop_c+1;
                fprintf('iteration %d of total %d\n', loop_c, it_total)
            end
        end

        function xf = bps_filter(obj,x) % Band pass filter
            fnyq   = obj.Empirical.fsamp/2;
            [b, a] = butter(obj.filter.order, [obj.filter.lp, obj.filter.hp]./fnyq);
            xf     = filtfilt(b,a,x);
        end


    end


end