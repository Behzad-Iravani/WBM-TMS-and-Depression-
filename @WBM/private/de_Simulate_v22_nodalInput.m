% -*- coding: 'UTF-8' -*-
function [X,Y,w] = de_Simulate_v22_nodalInput(obj)
% This function simulates time course using adaptive frequency model and
% heun method
% setup & initializing
%   Authors:
%           Neda Kaboodvand, n.kaboodvand@gmail.com 
%           Behzad Iravani, behzadiravani@gmail.com 
% This function is part of scripts for Macroscopic resting state model predicts
% theta burst stimulation response: a randomized trial


X    = single(zeros(size(obj.SC,1),obj.time.t_steps)); % pre-allocated signal of real component
Y    = single(zeros(size(obj.SC,1),obj.time.t_steps)); % pre-allocated signal of imaginary component
w    = single(zeros(size(obj.SC,1),obj.time.t_steps));
dsig = sqrt(obj.dt)*obj.random.beta; % precalculated timestep for noise

w(:,1) = obj.Empirical.freq'.*2*pi; %initializing with empirical BOLD freq
w(:,2) = w(:,1);

Diffx = X(:,1)-X(:,1)';
Diffy = Y(:,1)-Y(:,1)';

X(:,2) = X(:,1) + obj.dt*((obj.parameters.A-X(:,1).^2-Y(:,1).^2).*X(:,1) -w(:,2).*Y(:,1) +...
    obj.parameters.G.*sum(-1*obj.SC.*Diffx,2)) +dsig*obj.random.RR(:,1);
Y(:,2) = Y(:,1) + obj.dt*((obj.parameters.A-X(:,1).^2-Y(:,1).^2).*Y(:,1) +w(:,2).*X(:,1) +...
    obj.parameters.G.*sum(-1*obj.SC.*Diffy,2)) +dsig*obj.random.RR(:,1);

for t= 3:obj.time.t_steps
    % Heun integration
    % -- left slope
    [FX0,FY0,Fw0] = def_equ(obj, X(:,t-1),Y(:,t-1),w(:,t-1));
    w1 = w(:,t-1) + obj.dt*Fw0;
    X1 = X(:,t-1) + obj.dt*FX0  + dsig* obj.random.RR(:,t-1);
    Y1 = Y(:,t-1) + obj.dt*FY0  + dsig* obj.random.RR(:,t-1);
    % -- right slope
    [FX1,FY1,Fw1] = def_equ(obj, X1,Y1,w1);
    % new varaible
    w(:,t) = w(:,t-1) + obj.dt*(Fw0+Fw1)/2;
    X(:,t) = X(:,t-1) + obj.dt*(FX0+FX1)/2  + dsig* obj.random.RR(:,t-1);
    Y(:,t) = Y(:,t-1) + obj.dt*(FY0+FY1)/2  + dsig* obj.random.RR(:,t-1);
end

X          = X';
X(1:600,:) = [];

X = downsample(X ,fix((1/obj.Empirical.fsamp)/obj.dt));
% for c = 1:size(X,2)
% Xd(:,c)     = decimate(X(:,c) ,fix((1/obj.Empirical.fsamp)/obj.dt));%fix(1/obj.dt*obj.Empirical.fsamp)
% end
% X = Xd;
% clear Xd

Y          = Y';
Y(1:600,:) = [];
Y = downsample(Y ,fix((1/obj.Empirical.fsamp)/obj.dt));
% for c = 1:size(Y,2)
% Yd(:,c)      = decimate(Y(:,c),fix(1/obj.dt*obj.Empirical.fsamp));
% end
% Y = Yd;
% clear Yd

w          = w';
w(1:600,:) = [];
w = downsample(w ,fix((1/obj.Empirical.fsamp)/obj.dt));
% for c = 1:size(w,2)
% wd(:,c)     = decimate(w(:,c),fix(1/obj.dt*obj.Empirical.fsamp))';
% end
% w = wd;
% clear wd

    function [oX,oY,ow] = def_equ(obj, X,Y,w)

        %% Frequency
         
        
        ow     = (obj.parameters.F*w + ...
            obj.parameters.M.*(obj.SC...../repmat(sum(obj.SC),size(obj.SC,1),1))% normalize SC 
            *atan(Y./X))+ obj.Empirical.freq'.*2*pi);

        Diffx = X-X';
        Diffy = Y-Y';
        %% Amplitude
        oX = ((obj.parameters.A-X.^2-Y.^2).*X -...
            w.*Y + obj.parameters.G.*sum(-1*obj.SC...../repmat(sum(obj.SC),size(obj.SC,1),1))% normalize SC
            .*Diffx,2)) ;
        oY = ((obj.parameters.A-X.^2-Y.^2).*Y +...
            w.*X + obj.parameters.G.*sum(-1*obj.SC...../repmat(sum(obj.SC),size(obj.SC,1),1))% normalize SC
            .*Diffy,2));
    end
end

