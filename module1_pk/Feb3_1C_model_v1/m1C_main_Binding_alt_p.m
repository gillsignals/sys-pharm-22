clear all;
% Toy Model (i.e. generic model not specific to an actual drug):
%  specifically here, a one-compartment model with continuous infusion, 
%  binding to plasma protein, and different clearance rates for three molecules

% Here we show how to use a structure for parameters instead of passing a
% vector of values

%% Define parameters

% instead of defining each parameter as a separate variable...
% q = 1; % nmol/hr (drug input - continuous infusion)
% V = 1; % L (compartment volume)
% kcA = .1; % hr-1 (rate constant for elimination of molecule A)
% kcB = 0; % hr-1 (rate constant for elimination of molecule B)
% kcAB = .1; % hr-1 (rate constant for elimination of molecule AB)
% kon = .01; % nM-1 hr-1 (on-rate constant of A binding B)
% koff = .1; % hr-1 (off-rate constant of A dissociating from B)

% ...we can define a structure p of which each parameter is an element:
p.q = 1; % nmol/hr (drug input - continuous infusion)
p.V = 1; % L (compartment volume)
p.kcA = .1; % hr-1 (rate constant for elimination of molecule A)
p.kcB = 0; % hr-1 (rate constant for elimination of molecule B)
p.kcAB = .1; % hr-1 (rate constant for elimination of molecule AB)
p.kon = .01; % nM-1 hr-1 (on-rate constant of A binding B)
p.koff = .1; % hr-1 (off-rate constant of A dissociating from B)

% intial conditions (concentrations)
%   three molecules => three equations => three initial concentrations
%   order of elements is the same as order of equations
y0 = [0 100 0]'; % nM ([A B AB])

% now we don't need to assemble the parameters into a vector (but
% everything else is basically the same!)
% p = [q V kcA kcB kcAB kon koff]';

% assign some options values for the solver
options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);

%% RUN SIMULATIONS

% first simulation (based on parameter values set above)
% NOTE that p is still passed, but now is a structure instead of a vector;
% it's still being passed to the solver (ode45), which will in turn pass it
% to the equations function when it repeatedly calls that function
[T1,Y1] = ode45(@m1C_eqns_Binding_alt_p,[0 10],y0,options,p);
% calling the solver (in this case, ode45); to the solver, we pass
% the following arguments:
%   the name of the equations function (m1C_eqns_Binding_alt_p)
%   the time limits (0-10 hrs) 
%   initial conditions (y0)
%   ODE solver options (options)
%   parameter set (p) - note that this (and any other arguments you add after it)
%      gets passed through the solver into the equations function
% the ode45 function returns two things:
%   T1 = the vector of timesteps for which concentrations are returned to us
%   Y1 = the vector of concentrations at each timestep of T1
%        (since there are multiple equations, Y1 will have multiple columns)

% convert concentration to amount (nM -> nmoles)
TotalD1 = Y1*p.V ;

% second simulation - change values, run solver again
% p.q = 2;
% p.kcAB = .01;
p.kon = .1;

% because we updated the element of p directly, 
% we no longer need to re-assign values of p; p is already updated
% p = [q V kcA kcB kcAB kon koff]'; 

[T2,Y2] = ode45(@m1C_eqns_Binding_alt_p,[0 10],y0,options,p);
TotalD2 = Y2*p.V ;

% third simulation - change values, run solver again
% p.q = 3;
% p.kcAB = .001;
p.kon = 0.001;
[T3,Y3] = ode45(@m1C_eqns_Binding_alt_p,[0 10],y0,options,p);
TotalD3 = Y3*p.V ;

%% Visualize results

figure;
ax1=subplot(2,2,1);
plot(ax1,T1,Y1(:,1),'k',T1,Y1(:,3),'k-.',T2,Y2(:,1),'b',T2,Y2(:,3),'b-.',T3,Y3(:,1),'r',T3,Y3(:,3),'r-.','linewidth',3)
title(ax1,'Concentration of A(---),AB(- -)')
ylabel(ax1,'[D] (nM)')
xlabel(ax1,'time (hrs)')
% lgd = legend('1 (A)', '1 (AB)', '2 (A)',  '2 (AB)', '3 (A)', '3 (AB)');
% lgd.Location = 'northwest';
% lgd.Title.String = ['infusion rate' newline 'q (nmol/hr)'];
% lgd = legend('0.1 (A)', '0.1 (AB)', '0.01 (A)',  '0.01 (AB)', '0.001 (A)', '0.001 (AB)');
% lgd.Location = 'northwest';
% lgd.Title.String = ['AB clearance' newline 'rate constant' newline 'kcAB (hr^-^1)'];
lgd = legend('0.01 (A)', '0.01 (AB)', '0.1 (A)',  '0.1 (AB)', '0.001 (A)', '0.001 (AB)');
lgd.Location = 'northwest';
lgd.Title.String = ['A-B binding' newline 'rate constant' newline 'kon (nM^-^1 hr^-^1)'];

ax2=subplot(2,2,2);
plot(ax2,T1,Y1(:,2),'k',T2,Y2(:,2),'b',T3,Y3(:,2),'r','linewidth',3)
title(ax2,'Concentration of B')
ylabel(ax2,'[D] (nM)')
xlabel(ax2,'time (hrs)')
% lgd = legend('1', '2', '3');
% lgd.Location = 'southwest';
% lgd.Title.String = ['infusion rate' newline 'q (nmol/hr)'];
% lgd = legend('.1', '0.01', '0.001');
% lgd.Location = 'southwest';
% lgd.Title.String = ['AB clearance' newline 'rate constant' newline 'kcAB (hr^-^1)'];
lgd = legend('0.01', '0.1', '0.001');
lgd.Location = 'southwest';
lgd.Title.String = ['A-B binding' newline 'rate constant' newline 'kon (nM^-^1 hr^-^1)'];

ax3=subplot(2,2,3);
plot(ax3,T1,TotalD1(:,1),'k',T2,TotalD2(:,1),'b',T3,TotalD3(:,1),'r','linewidth',3)
title(ax3,'Total Amount of Free Drug in Compartment')
ylabel(ax3,'Total Drug (nmol)')
xlabel(ax3,'time (hrs)')
% lgd = legend('1', '2', '3');
% lgd.Location = 'northwest';
% lgd.Title.String = ['infusion rate' newline 'q (nmol/hr)'];
% lgd = legend('.1', '0.01', '0.001');
% lgd.Location = 'northwest';
% lgd.Title.String = ['AB clearance' newline 'rate constant' newline 'kcAB (hr^-^1)'];
lgd = legend('0.01', '0.1', '0.001');
lgd.Location = 'northwest';
lgd.Title.String = ['A-B binding' newline 'rate constant' newline 'kon (nM^-^1 hr^-^1)'];

ax4=subplot(2,2,4);
plot(ax4,T1,TotalD1(:,1)+TotalD1(:,3),'k',T2,TotalD2(:,1)+TotalD2(:,3),'b',T3,TotalD3(:,1)+TotalD3(:,3),'r','linewidth',3)
title(ax4,'Total Amount of Total Drug in Compartment')
ylabel(ax4,'Total Drug (nmol)')
xlabel(ax4,'time (hrs)')
% lgd = legend('1', '2', '3');
% lgd.Location = 'northwest';
% lgd.Title.String = ['infusion rate' newline 'q (nmol/hr)'];
% lgd = legend('.1', '0.01', '0.001');
% lgd.Location = 'northwest';
% lgd.Title.String = ['AB clearance' newline 'rate constant' newline 'kcAB (hr^-^1)'];
lgd = legend('0.01', '0.1', '0.001');
lgd.Location = 'northwest';
lgd.Title.String = ['A-B binding' newline 'rate constant' newline 'kon (nM^-^1 hr^-^1)'];

% Note - we will explore more compact methods of writing code to run 
% multiple simulations, and more compact methods of writing code to
% visualize those simulations; we've written these out completely to give
% you a clearer understanding of the code before we get into being more
% efficient with number of lines of code!