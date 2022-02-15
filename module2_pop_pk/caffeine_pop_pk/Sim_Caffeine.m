function [out,out2] = Sim_Caffeine(w,x,y,z,f,visFlag,TimeLen); 
% NOTE - for one of the cases (case 5), we need two outputs, for the others
% we only need one. So we can set up the function with two, and just return
% a dummy value (0) for the second one. The outputs in each case are different too,
% so the case flag ('f') will determine what gets returned below.

% This function runs the caffeine/paraxanthine simulations, according to
% equations described in CaffParax_eqns.m

%% PARAMETERS
% Model parameters for array p
q = 0; % mg/hr
V1 = x; % L

halflife = y; % hr
kc1 = 0.693/halflife; % hr-1
kc2 = kc1/2; % hr-1 paraxanthine clearance slower

kcp = 0;  % hr-1 conversion rate of caffeine to paraxanthine

abshalflife = w; % hr
ka = 0.693/abshalflife; % hr-1

% Initial concentrations - caffeine in gut
D0 = z; % mg

% TimeLen = 2; %hrs

y0 = [0 0 0 0 D0]'; % nM - Initial conditions vector
  % 1 = caffeine in body (mg/L)
  % 2 = paraxanthine in body (mg/L)
  % 3 = caffeine in degr (mg)
  % 4 = paraxanthine in degr (mg)
  % 5 = caffeine in gut (mg)
  
p = [q V1 kc1 kc2 kcp ka]'; % parameter array

%% RUN SIMULATION

options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
[T1,Y1] = ode45(@CaffParax_eqns,[0 TimeLen],y0,options,p);

%% CALCULATE OUTPUTS

DrugInBolus = ones(length(T1),1)*(D0) ; % this is in mg
TotalFreeD(:,1) = Y1(:,1)*V1; % this is in mg (mg/L * L)
TotalFreeD(:,2) = Y1(:,2)*V1; % this is in mg (mg/L * L)
DrugIn = q*T1 + DrugInBolus ; % cumulative drug into system - both in mg
DrugOut = Y1(:,3) +  Y1(:,4) ; % cumulative drug eliminated from system - this is in mg already
BalanceD1 = DrugIn - DrugOut - TotalFreeD(:,1) - TotalFreeD(:,2) ; %(zero = balance)

% Outputs depend on which case is being run

switch f
    case {1,2,6}
        % OUTPUT: peak concentration of caffeine (single number)
        out = max(Y1(:,1));
        out2 = 0; % dummy
    case {3,4}
        % OUTPUT: concentrations of caffeine over time (array)
        out = Y1(:,1);
        out2 = 0; % dummy
    case {5,12,13}
        % OUTPUT: list of timepoints and concentrations of caffeine over time
        out = T1(:);
        out2 = Y1(:,1);
    case {11}
        % OUTPUT: list of timepoints and concentrations of caffeine over time
        out = T1(:);
        out2 = Y1(:,5);
end

%% VISUALIZATION

figure('visible',visFlag); % only show figures if flag is 'on'

ax1=subplot(3,2,1);
plot(ax1,T1,Y1(:,1),'k','linewidth',3)
title(ax1,'Concentration of Caffeine in Blood')
ylabel(ax1,'[Caff] (mg/L)')
xlabel(ax1,'time (hrs)')

ax2=subplot(3,2,2);
plot(ax2,T1,Y1(:,2),'k','linewidth',3)
title(ax2,'Concentration of Paraxanthine in Blood')
ylabel(ax2,'[Parax] (mg/L)')
xlabel(ax2,'time (hrs)')

ax3=subplot(3,2,3);
plot(ax3,T1,Y1(:,5),'k','linewidth',3)
title(ax3,'Amount of Caffeine in Gut')
ylabel(ax3,'Caffeine (mg)')
xlabel(ax3,'time (hrs)')

ax4=subplot(3,2,5);
plot(ax4,T1,Y1(:,3),'k','linewidth',3)
title(ax4,'Amount of Cleared Caffeine')
ylabel(ax4,'Caffeine (mg)')
xlabel(ax4,'time (hrs)')

ax5=subplot(3,2,6);
plot(ax5,T1,Y1(:,4),'k','linewidth',3)
title(ax5,'Amount of Cleared Paraxanthine')
ylabel(ax5,'Paraxanthine (mg)')
xlabel(ax5,'time (hrs)')
