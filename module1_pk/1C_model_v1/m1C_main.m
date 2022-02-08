clear all;
% Toy Model (i.e. generic model not specific to an actual drug):
%  specifically here, a one-compartment model with infusion and clearance

%% DEFINE PARAMETER VALUES
q = 1; % nmol/hr (drug input - continuous infusion)
V = 1; % L (compartment volume)
k = 1; % hr-1 (rate constant for elimination)

% intial conditions (concentrations); y0 is a vector if multiple equations 
y0 = 0; % nM 

% Assemble parameters into a vector, to simplify passing the values to the solver
% note that the order of the parameters here needs to be consistent with the order in the eqns file
p = [q V k]'; 

% assign some options values for the solver
options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);


%% RUN SIMULATIONS
% Note that running the simulations happens in just one line of code!

% first simulation
[T1,Y1] = ode45(@m1C_eqns,[0 10],y0,options,p);
% calling the solver (in this case, ode45); to the solver, we pass
% the following arguments:
%   the name of the equations function (m1C_eqns)
%   the time limits (0-10 hrs) 
%   initial conditions (y0)
%   ODE solver options (options)
%   parameter array (p) - note that this (and any other arguments you add after it)
%      gets passed through the solver into the equations function
% the ode45 function returns two things:
%   T1 = the vector of timesteps for which concentrations are returned to us
%   Y1 = the vector of concentrations at each timestep of T1
%        (if there are multiple equations, Y1 will have multiple columns)

% convert concentration to amount (nM -> nmoles)
TotalD1 = Y1*V ;  % concentration * volume = amount

% second simulation - change parameter values, run solver again
q = 2;
p = [q V k]'; % the vector p doesn't automatically update, need to reassign updated values to it
[T2,Y2] = ode45(@m1C_eqns,[0 10],y0,options,p);
TotalD2 = Y2*V ;

% third simulation - change parameter values, run solver again
q = 3;
p = [q V k]'; % the vector p doesn't automatically update, need to reassign updated values to it
[T3,Y3] = ode45(@m1C_eqns,[0 10],y0,options,p);
TotalD3 = Y3*V ;

%% VISUALIZE RESULTS

figure;
ax1=subplot(1,2,1);
plot(ax1,T1,Y1(:,1),'k',T2,Y2(:,1),'b',T3,Y3(:,1),'r','linewidth',3)
title(ax1,'Concentration of Drug in Compartment')
ylabel(ax1,'[D] (nM)')
xlabel(ax1,'time (hrs)')
lgd = legend('1', '2', '3');
lgd.Location = 'northeast';
lgd.Title.String = ['infusion rate' newline 'q (nmol/hr)'];
% lgd.Title.String = ['Volume' newline 'V (L)'];
% lgd.Title.String = ['clearance' newline 'rate constant' newline 'k (1/hr)'];

ax2=subplot(1,2,2);
plot(ax2,T1,TotalD1(:,1),'k',T2,TotalD2(:,1),'b',T3,TotalD3(:,1),'r','linewidth',3)
title(ax2,'Total Amount of Drug in Compartment')
ylabel(ax2,'Total Drug (nmol)')
xlabel(ax2,'time (hrs)')
lgd = legend('1', '2', '3');
lgd.Location = 'northeast';
lgd.Title.String = ['infusion rate' newline 'q (nmol/hr)'];
% lgd.Title.String = ['Volume' newline 'V (L)'];
% lgd.Title.String = ['clearance' newline 'rate constant' newline 'k (1/hr)'];
