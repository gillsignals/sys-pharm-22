clear all;
% Toy Model (i.e. generic model not specific to an actual drug):
%  specifically here, a one-compartment model with two molecules
%  (a drug or prodrug, and its metabolite) and multiple processes
%  (infusion, clearance, and metabolism)

%% DEFINE PARAMETER VALUES
q = 1; % nmol/hr (drug input - continuous infusion)
V = 1; % L (compartment volume)
kc1 = .1; % hr-1 (rate constant for elimination of molecule A)
kc2 = .1; % hr-1 (rate constant for elimination of molecule B)
kab = 1; % hr-1 (rate constant for first-order conversion of A => B)

% intial conditions (concentrations); y0 is a vector here 
%   (two molecules => two equations => two initial concentrations)
y0 = [0 0]'; % nM

% Assemble parameters into a vector, to simplify passing the values to the solver
% note that the order of the parameters here needs to be consistent with the order in the eqns file
p = [q V kc1 kc2 kab]'; 

% assign some options values for the solver
options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);

%% RUN SIMULATIONS
% Note that running the simulations happens in just one line of code!

% first simulation
[T1,Y1] = ode45(@m1C_eqns_Metab,[0 10],y0,options,p);
% calling the solver (in this case, ode45); to the solver, we pass
% the following arguments:
%   the name of the equations function (m1C_eqns_Metab)
%   the time limits (0-10 hrs) 
%   initial conditions (y)
%   ODE solver options (options)
%   parameter array (p) - note that this (and any other arguments you add after it)
%      gets passed through the solver into the equations function
% the ode45 function returns two things:
%   T1 = the vector of timesteps for which concentrations are returned to us
%   Y1 = the vector of concentrations at each timestep of T1
%        (since there are multiple equations, Y1 will have multiple columns)

% convert concentration to amount (nM -> nmoles)
TotalD1 = Y1*V ; % concentration * volume = amount

% second simulation - change parameter values, run solver again
% q = 2;
% kc2 = .01;
kab = .1; 
p = [q V kc1 kc2 kab]'; % the vector p doesn't automatically update, need to reassign updated values to it
[T2,Y2] = ode45(@m1C_eqns_Metab,[0 10],y0,options,p);
TotalD2 = Y2*V ;

% third simulation - change parameter values, run solver again
% q = 3;
% kc2 = .001;
kab = .01; 
p = [q V kc1 kc2 kab]'; % the vector p doesn't automatically update, need to reassign updated values to it
[T3,Y3] = ode45(@m1C_eqns_Metab,[0 10],y0,options,p);
TotalD3 = Y3*V ;

%% VISUALIZE RESULTS

figure;
ax1=subplot(2,2,1);
plot(ax1,T1,Y1(:,1),'k',T2,Y2(:,1),'b',T3,Y3(:,1),'r','linewidth',3)
title(ax1,'Concentration of A')
ylabel(ax1,'[D] (nM)')
xlabel(ax1,'time (hrs)')
% lgd = legend('1', '2', '3');
% lgd.Location = 'northwest';
% lgd.Title.String = ['infusion rate' newline 'q (nmol/hr)'];
% lgd = legend('.1', '0.01', '0.001');
% lgd.Location = 'northwest';
% lgd.Title.String = ['B clearance' newline 'rate constant' newline 'kc2 (hr^-^1)'];
lgd = legend('1', '0.1', '0.01');
lgd.Location = 'northwest';
lgd.Title.String = ['A->B conversion' newline 'rate constant' newline 'kab (hr^-^1)'];

ax2=subplot(2,2,2);
plot(ax2,T1,Y1(:,2),'k',T2,Y2(:,2),'b',T3,Y3(:,2),'r','linewidth',3)
title(ax2,'Concentration of B')
ylabel(ax2,'[D] (nM)')
xlabel(ax2,'time (hrs)')
% lgd = legend('1', '2', '3');
% lgd.Location = 'northwest';
% lgd.Title.String = ['infusion rate' newline 'q (nmol/hr)'];
% lgd = legend('.1', '0.01', '0.001');
% lgd.Location = 'northwest';
% lgd.Title.String = ['B clearance' newline 'rate constant' newline 'kc2 (hr^-^1)'];
lgd = legend('1', '0.1', '0.01');
lgd.Location = 'northwest';
lgd.Title.String = ['A->B conversion' newline 'rate constant' newline 'kab (hr^-^1)'];

ax3=subplot(2,2,3);
plot(ax3,T1,TotalD1(:,1),'k',T2,TotalD2(:,1),'b',T3,TotalD3(:,1),'r','linewidth',3)
title(ax3,'Total Amount of Free Drug (A) in Compartment')
ylabel(ax3,'Total Drug (nmol)')
xlabel(ax3,'time (hrs)')
% lgd = legend('1', '2', '3');
% lgd.Location = 'northwest';
% lgd.Title.String = ['infusion rate' newline 'q (nmol/hr)'];
% lgd = legend('.1', '0.01', '0.001');
% lgd.Location = 'northwest';
% lgd.Title.String = ['B clearance' newline 'rate constant' newline 'kc2 (hr^-^1)'];
lgd = legend('1', '0.1', '0.01');
lgd.Location = 'northwest';
lgd.Title.String = ['A->B conversion' newline 'rate constant' newline 'kab (hr^-^1)'];

ax4=subplot(2,2,4);
plot(ax4,T1,TotalD1(:,1)+TotalD1(:,2),'k',T2,TotalD2(:,1)+TotalD2(:,2),'b',T3,TotalD3(:,1)+TotalD3(:,2),'r','linewidth',3)
title(ax4,'Total Amount of Total Drug (A+B) in Compartment')
ylabel(ax4,'Total Drug (nmol)')
xlabel(ax4,'time (hrs)')
% lgd = legend('1', '2', '3');
% lgd.Location = 'northwest';
% lgd.Title.String = ['infusion rate' newline 'q (nmol/hr)'];
% lgd = legend('.1', '0.01', '0.001');
% lgd.Location = 'northwest';
% lgd.Title.String = ['B clearance' newline 'rate constant' newline 'kc2 (hr^-^1)'];
lgd = legend('1', '0.1', '0.01');
lgd.Location = 'northwest';
lgd.Title.String = ['A->B conversion' newline 'rate constant' newline 'kab (hr^-^1)'];

% Note - we will explore more compact methods of writing code to run 
% multiple simulations, and more compact methods of writing code to
% visualize those simulations; we've written these out completely to give
% you a clearer understanding of the code before we get into being more
% efficient with number of lines of code!
