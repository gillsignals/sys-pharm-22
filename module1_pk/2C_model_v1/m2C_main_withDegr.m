close all;
clear all;
% Toy Model (i.e. generic model not specific to an actual drug):
%  specifically here, a two-compartment model but with an additional 3rd
%  'virtual compartment' to keep track of drug cleared from the system. 

%% DEFINE PARAMETER VALUES
p.q = 1; % nmol/hr (drug input - continuous infusion)
p.V1 = 1; % L (central compartment volume)
p.V2 = 1; % L (peripheral compartment volume)
p.kc1 = 1; % hr-1 (rate constant for elimination from central compartment)
p.kc2 = 0; % hr-1 (rate constant for elimination from peripheral compartment)
p.k12 = 1; % hr-1 (rate constant for transport from central to peripheral compartment)
p.k21 = 1; % hr-1 (rate constant for transport from peripheral to central compartment)

% intial conditions (concentrations and/or amounts)
%  y0 is a vector because there are three quantities 
%  (two concentrations and one amount) being simulated
y0 = [0 0 0]'; 
% 1st element of y0 = drug in central compartment; unit is nM
% 2nd element of y0 = drug in peripheral compartment; unit is nM
% 3rd element of y0 = drug in cleared compartment; unit is nmol (amount)

% assign some options values for the solver
options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);

%% RUN SIMULATIONS

% first simulation (based on parameter values defined above)
[T1,Y1] = ode45(@m2C_eqns_withDegr,[0:(1/60):10],y0,options,p);
% calling the solver (in this case, ode45); to the solver, we pass
% the following arguments:
%   the name of the equations function (m2C_eqns_withDegr)
%   the time limits (0-10 hrs) 
%   initial conditions (y0)
%   ODE solver options (options)
%   parameter set (p) - note that this (and any other arguments you add after it)
%      gets passed through the solver into the equations function
% the ode45 function returns two things:
%   T1 = the vector of timesteps for which concentrations are returned to us
%   Y1 = the vector of concentrations or amounts at each timestep of T1
%        (since there are three equations, Y1 will have three columns)

%% CALCULATE MASS BALANCE COMPONENTS
CurrentDrug(:,1) = Y1(:,1)*p.V1; % Total drug in central compartment at time t
CurrentDrug(:,2) = Y1(:,2)*p.V2; % Total drug in peripheral compartment at time t
InitialDrug = y0(1)*p.V1 + y0(2)*p.V2 ; % Total drug in both compartments at time zero 
DrugIn = p.q*T1 ;   % Cumulative drug into system (assumes constant q)
DrugOut = Y1(:,3) ; % Cumulative drug eliminated from system
% Calculate mass balance error (zero = balance)
BalanceD = DrugIn - DrugOut - CurrentDrug(:,1) - CurrentDrug(:,2) + InitialDrug ; %(zero = balance)

%p.q = 2;
%[T2,Y2] = ode45(@m2C_eqns_withDegr,[0 10],y0,options,p);

%p.q = 3;
%[T3,Y3] = ode45(@m2C_eqns_withDegr,[0 10],y0,options,p);

%% VISUALIZE RESULTS

fig1 = figure;
ax1=subplot(2,2,1);
plot(ax1,T1,Y1(:,1),'k',T1,Y1(:,2),'r','linewidth',3)
title(ax1,'Concentration of Drug in Compartments')
ylabel(ax1,'[D] (nM)')
xlabel(ax1,'time (hrs)')
lgd = legend('central', 'peripheral');
lgd.Location = 'southeast';
lgd.Title.String = ['Compartment'];

ax2=subplot(2,2,2);
plot(ax2,T1,CurrentDrug(:,1),'k',T1,CurrentDrug(:,2),'r','linewidth',3)
title(ax2,'Total Amount of Drug in Compartments')
ylabel(ax2,'Total Drug (nmol)')
xlabel(ax2,'time (hrs)')
lgd = legend('central', 'peripheral');
lgd.Location = 'southeast';
lgd.Title.String = ['Compartment'];

ax3=subplot(2,2,3);
plot(ax3,T1,DrugIn,'b-',T1,DrugOut,'b-.','linewidth',3)
title(ax3,'Cumulative drug in or out of system') 
ylabel(ax3,['Cumulative Drug' newline 'in/out of system (nmol)'])
xlabel(ax3,'time (hrs)')
lgd = legend('Input', 'Output');
lgd.Location = 'southeast'; 
lgd.Title.String = ['Drug (nmol)'];

ax4=subplot(2,2,4);
plot(ax4,T1,BalanceD,'m','linewidth',3)
title(ax4,'Molecular Balance') %(zero = balance)
ylabel(ax4,'Balance of Drug (nmol)')
xlabel(ax4,'time (hrs)')

%% EXPORT VISUALIZATION
set(fig1,'Position',[0 0 600 450])
exportgraphics(fig1, "Fig1_TwoComptModel.png",'Resolution',300); 
