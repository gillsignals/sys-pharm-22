clear all;
% Toy one-compartment model 

%% Define parameters
q = 1; % nmol/hr (drug input - continuous infusion)
V = 1; % L (compartment volume)
k = 1; % hr-1 (rate constant for elimination)
y0 = 0; % nM (initial concentration; is a vector if multiple equations)
p = [q V k]'; % (parameter vector for equations - passed thru solver)

options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
% assigning some options values for the solver

%% Run simulations

% first simulation
[T1,Y1] = ode45(@m1C_eqns,[0 10],y0,options,p);
% calling the solver (in this case, ode45); pass the equations function,
% the time limits (0-10 hrs), initial conditions, ODE options, and the
% parameter array p

TotalD1 = Y1*V ;
% convert concentration to amount (nM -> nmoles)

% second simulation - change values, run solver again
q = 2;
p = [q V k]'; % p doesn't automatically update, need to reassign it
[T2,Y2] = ode45(@m1C_eqns,[0 10],y0,options,p);
TotalD2 = Y2*V ;

% third simulation - change values, run solver again
q = 3;
p = [q V k]';
[T3,Y3] = ode45(@m1C_eqns,[0 10],y0,options,p);
TotalD3 = Y3*V ;

%% Plot results

figure;
ax1=subplot(1,2,1);
plot(ax1,T1,Y1(:,1),'k',T2,Y2(:,1),'b',T3,Y3(:,1),'r','linewidth',3)
title(ax1,'Concentration of Drug in Compartment')
ylabel(ax1,'[D] (nM)')
xlabel(ax1,'time (hrs)')

ax2=subplot(1,2,2);
plot(ax2,T1,TotalD1(:,1),'k',T2,TotalD2(:,1),'b',T3,TotalD3(:,1),'r','linewidth',3)
title(ax2,'Total Amount of Drug in Compartment')
ylabel(ax2,'Total Drug (nmol)')
xlabel(ax2,'time (hrs)')
