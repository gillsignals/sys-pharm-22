close all;
clear all;
% Toy Model (i.e. generic model not specific to an actual drug):
%  specifically here, a one-compartment model but with an additional 
%  'virtual compartment' to keep track of drug cleared from the system. 
% We can approach this virtual compartment in (at least) two different ways: 

% (1) The virtual compartment tracks the cleared material as a
% CONCENTRATION (e.g. ng/ml, nM); in this case, the virtual clearance 
% compartment needs to have a volume associated with it; 
% we typically assume (for simplicity) that the cleared compartment 
% has the same volume as the central compartment. If the volume is
% different, then we would need to include volume corrections. Since the
% volume is arbitrary, simple to make it easier on ourselves in this way.

% (2) The virtual compartment tracks the cleared material as an AMOUNT
% (e.g. ng, nmol); in this case, we don't need a volume associated with the
% compartment, but we do need to convert the rate of clearance from
% concentration/time to amount/time; and we need to be careful in the mass
% balance because some of our variables have concentration units and some
% have amount units.

% On balance, approach (2) is probably more straightforward, since it
% doesn't lead to thinking about the virtual compartment as any sort of
% physical space; it's just a virtual bucket for collecting the cleared
% material.

%% DEFINE PARAMETER VALUES
p.q = 1; % nmol/hr (drug input - continuous infusion)
p.V = 1; % L (compartment volume)
p.k = 0.1; % hr-1 (rate constant for elimination)

% Select method for virtual compartment
p.cleared = 'conc';
% p.cleared = 'amnt';
%
% These are the two potential approaches for dealing with the virtual
% compartment. This is included in p so that the information about which 
% approach to use is passed to the equations file. You should only need to
% change the value of p.cleared, and later code will take care of the rest.
% Examine how the code does this, because this is a useful format you will
% be able to use later in the course.

% Select parameter to simulate variation
paramToVary = 'q';
% paramToVary = 'V';
% paramToVary = 'k';
% paramToVary = 'D0';
% 
% We can make our code slightly more flexible if we don't have to keep
% commenting and uncommenting different lines - for example, to change from
% varying q to varying k. So here we only have to uncomment one of the four
% "paramToVary" lines above and the code will take care of the rest below
% (including putting the right legend on the graphs). We can be even more
% efficient than this (for example by putting the visualizations into
% functions) but this will be an improvement. You should read the code
% carefully and learn how to apply this kind of thinking in your own codes.


D0 = 0; % nmol (bolus dose, if there is one)

% intial conditions (concentrations and/or amounts)
%  y0 is a vector because there are two concentrations or amounts tracked
y0 = [D0/p.V 0]'; 
% 1st element of y0 = drug in central compartment; unit is nM
% 2nd element of y0 = drug in cleared compartment; unit depends on 
%   choice of approach - nM for concentration, nmol for amount.

% assign some options values for the solver
options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);

%% RUN SIMULATIONS

% first simulation (based on parameter values defined above)
[T1,Y1] = ode45(@m1C_eqns_withDegr,[0:(1/60):10],y0,options,p);
% calling the solver (in this case, ode45); to the solver, we pass
% the following arguments:
%   the name of the equations function (m1C_eqns_withDegr)
%   the time limits (0-10 hrs) 
%   initial conditions (y0)
%   ODE solver options (options)
%   parameter set (p) - note that this (and any other arguments you add after it)
%      gets passed through the solver into the equations function
% the ode45 function returns two things:
%   T1 = the vector of timesteps for which concentrations are returned to us
%   Y1 = the vector of concentrations or amounts at each timestep of T1
%        (since there are two equations, Y1 will have two columns)

%% CALCULATE MASS BALANCE COMPONENTS
InitialDrug1 =   y0(1)*p.V ; % Total drug in central compartment at time zero 
CurrentDrug1 = Y1(:,1)*p.V ; % Total drug in central compartment at time t
% NOTE on syntax:
% Y1(:,1) means the all the values in the first column of Y1; the colon
% includes all values (https://www.mathworks.com/help/matlab/ref/colon.html)

DrugIn1 = p.q*T1 ; % Cumulative drug into system (assumes constant q)

if p.cleared == 'conc'
    DrugOut1 = Y1(:,2)*p.V ; % Cumulative drug eliminated from system
elseif p.cleared == 'amnt'
    DrugOut1 = Y1(:,2)     ; % Cumulative drug eliminated from system
else
    disp('error - no clearance compartment approach defined')
    return
end
% NOTE how we want amount for the mass balance, so we need to make sure to
% multiply by volume if using concentration in the virtual compartment

% Calculate mass balance error (zero = balance)
BalanceD1 = DrugIn1 - DrugOut1 - CurrentDrug1 + InitialDrug1; 
% Calculate relative mass balance error (zero = balance)
%   (based on initial drug and total drug added)
BalanceD1norm = BalanceD1./(DrugIn1+InitialDrug1) ; 

% second simulation - change parameter values, run solver again
switch paramToVary 
    case "q"
        p.q = 2;
    case "V"
        p.V = 2; 
    case "k"
        p.k = .01;
    case "D0"
        D0 = 1; 
end
y0 = [D0/p.V 0]'; % initial conditions
[T2,Y2] = ode45(@m1C_eqns_withDegr,[0:(1/60):10],y0,options,p);

% second simulation - calculate mass balance again
InitialDrug2 =   y0(1)*p.V ; % Total drug in central compartment at time zero 
CurrentDrug2 = Y2(:,1)*p.V ; % Total drug in central compartment at time t
DrugIn2 = p.q*T2 ; % Cumulative drug into system (assumes constant q)
if p.cleared == 'conc'
    DrugOut2 = Y2(:,2)*p.V ; % Cumulative drug eliminated from system
elseif p.cleared == 'amnt'
    DrugOut2 = Y2(:,2)     ; % Cumulative drug eliminated from system
else
    disp('error - no clearance compartment approach defined')
end

BalanceD2 = DrugIn2 - DrugOut2 - CurrentDrug2 + InitialDrug2; 
BalanceD2norm = BalanceD2./(DrugIn2+InitialDrug2) ; %(zero = balance)

% third simulation - change parameter values, run solver again
switch paramToVary 
    case "q"
        p.q = 3;
    case "V"
        p.V = 3; 
    case "k"
        p.k = .001;
    case "D0"
        D0 = 2; 
end
y0 = [D0/p.V 0]'; % nM
[T3,Y3] = ode45(@m1C_eqns_withDegr,[0:(1/60):10],y0,options,p);

% third simulation - calculate mass balance again
InitialDrug3 =   y0(1)*p.V ; % Total drug in central compartment at time zero 
CurrentDrug3 = Y3(:,1)*p.V ; % Total drug in central compartment at time t
DrugIn3 = p.q*T3 ; % Cumulative drug into system (assumes constant q)
if p.cleared == 'conc'
    DrugOut3 = Y3(:,2)*p.V ; % Cumulative drug eliminated from system
elseif p.cleared == 'amnt'
    DrugOut3 = Y3(:,2)     ; % Cumulative drug eliminated from system
else
    disp('error - no clearance compartment approach defined')
end

BalanceD3 = DrugIn3 - DrugOut3 - CurrentDrug3 + InitialDrug3 ; 
BalanceD3norm = BalanceD3./(DrugIn3+InitialDrug3) ; %(zero = balance)

%% VISUALIZE RESULTS

fig1 = figure;
ax1=subplot(3,2,1);
plot(ax1,T1,Y1(:,1),'k',T2,Y2(:,1),'b',T3,Y3(:,1),'r','linewidth',3)
title(ax1,'Concentration of Drug in Central Compartment')
ylabel(ax1,'[D] (nM)')
xlabel(ax1,'time (hrs)')
switch paramToVary
    case "q"
        lgd = legend('1', '2', '3');
        lgd.Location = 'northwest';
        lgd.Title.String = ['infusion rate' newline 'q (nmol/hr)'];
    case "V"
        lgd = legend('1', '2', '3');
        lgd.Location = 'northwest';
        lgd.Title.String = ['Volume' newline 'V (L)'];
    case "k"
        lgd = legend('0.1', '0.01', '0.001');
        lgd.Location = 'northwest';
        lgd.Title.String = ['clearance' newline 'rate constant' newline 'k (1/hr)'];
    case "D0"
        lgd = legend('0', '1', '2');
        lgd.Location = 'northwest';
        lgd.Title.String = ['Drug Bolus' newline 'D0 (nM)'];
end

ax2=subplot(3,2,2);
plot(ax2,T1,CurrentDrug1,'k',T2,CurrentDrug2,'b',T3,CurrentDrug3,'r','linewidth',3)
title(ax2,'Total Amount of Drug in Central Compartment')
ylabel(ax2,'Total Drug (nmol)')
xlabel(ax2,'time (hrs)')
switch paramToVary
    case "q"
        lgd = legend('1', '2', '3');
        lgd.Location = 'northwest';
        lgd.Title.String = ['infusion rate' newline 'q (nmol/hr)'];
    case "V"
        lgd = legend('1', '2', '3');
        lgd.Location = 'northwest';
        lgd.Title.String = ['Volume' newline 'V (L)'];
    case "k"
        lgd = legend('0.1', '0.01', '0.001');
        lgd.Location = 'northwest';
        lgd.Title.String = ['clearance' newline 'rate constant' newline 'k (1/hr)'];
    case "D0"
        lgd = legend('0', '1', '2');
        lgd.Location = 'northwest';
        lgd.Title.String = ['Drug Bolus' newline 'D0 (nM)'];
end

ax3=subplot(3,2,3);
plot(ax3,T1,DrugIn1,'k-',T2,DrugIn2,'b-',T3,DrugIn3,'r-',T1,DrugOut1,'k-.',T2,DrugOut2,'b-.',T3,DrugOut3,'r-.','linewidth',3)
title(ax3,'Cumulative Input (---)/Output(- -)') 
ylabel(ax3,['Cumulative Drug' newline 'in/out of system (nmol)'])
xlabel(ax3,'time (hrs)')
switch paramToVary
    case "q"
        lgd = legend('1', '2', '3');
        lgd.Location = 'northwest';
        lgd.Title.String = ['infusion rate' newline 'q (nmol/hr)'];
    case "V"
        lgd = legend('1', '2', '3');
        lgd.Location = 'northwest';
        lgd.Title.String = ['Volume' newline 'V (L)'];
    case "k"
        lgd = legend('0.1', '0.01', '0.001');
        lgd.Location = 'northwest';
        lgd.Title.String = ['clearance' newline 'rate constant' newline 'k (1/hr)'];
    case "D0"
        lgd = legend('0', '1', '2');
        lgd.Location = 'northwest';
        lgd.Title.String = ['Drug Bolus' newline 'D0 (nM)'];
end

ax4=subplot(3,2,4);
plot(ax4,T1,BalanceD1,'k-',T2,BalanceD2,'r-',T3,BalanceD3,'b-','linewidth',3)
title(ax4,'Molecular Balance') %(zero = balance)
ylabel(ax4,'Balance of Drug (nmol)')
xlabel(ax4,'time (hrs)')
switch paramToVary
    case "q"
        lgd = legend('1', '2', '3');
        lgd.Location = 'northwest';
        lgd.Title.String = ['infusion rate' newline 'q (nmol/hr)'];
    case "V"
        lgd = legend('1', '2', '3');
        lgd.Location = 'northwest';
        lgd.Title.String = ['Volume' newline 'V (L)'];
    case "k"
        lgd = legend('0.1', '0.01', '0.001');
        lgd.Location = 'northwest';
        lgd.Title.String = ['clearance' newline 'rate constant' newline 'k (1/hr)'];
    case "D0"
        lgd = legend('0', '1', '2');
        lgd.Location = 'northwest';
        lgd.Title.String = ['Drug Bolus' newline 'D0 (nM)'];
end

ax5=subplot(3,2,[5 6]); % "[5 6]" notation stretches the subplot across two panels
plot(ax5,T1,BalanceD1norm,'k-',T2,BalanceD2norm,'r-',T3,BalanceD3norm,'b-','linewidth',3)
title(ax5,'Molecular Balance') %(zero = balance)
ylabel(ax5,['Balance of Drug' newline ' (frac of drug in)'])
xlabel(ax5,'time (hrs)')
switch paramToVary
    case "q"
        lgd = legend('1', '2', '3');
        lgd.Location = 'northwest';
        lgd.Title.String = ['infusion rate' newline 'q (nmol/hr)'];
    case "V"
        lgd = legend('1', '2', '3');
        lgd.Location = 'northwest';
        lgd.Title.String = ['Volume' newline 'V (L)'];
    case "k"
        lgd = legend('0.1', '0.01', '0.001');
        lgd.Location = 'northwest';
        lgd.Title.String = ['clearance' newline 'rate constant' newline 'k (1/hr)'];
    case "D0"
        lgd = legend('0', '1', '2');
        lgd.Location = 'northwest';
        lgd.Title.String = ['Drug Bolus' newline 'D0 (nM)'];
end

%% EXPORT VISUALIZATION
% Rather than copy figures from Matlab's windows, we can exert more control
% over the size and resolution of the figure(s) and export them to files
% (named appropriately and of a file format we choose)
set(fig1,'Position',[0 0 600 600])
exportgraphics(fig1, strcat("Fig1_varying_",paramToVary,".png"),'Resolution',300); 
% Notice how the filename is different depending on which parameter is
% being varied across the three simulations