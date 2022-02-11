clear all;
% Two-compartment model of acetaminophen. This code runs two simulations -
% first, using parameters appropriate for the softgel formulation; 
% second, using parameters appropriate for the IV delivery.
% Then, it compares the AUC values for the two delivery modes to estimate
% bioavailability for the oral delivery

%% PARAMETER VALUES - SOFTGEL
p.q = 0; % ug/hr  (infusion rate to central compartment; for IV delivery)
p.V1 = 5000; % mL (Central compartment volume)
p.V2 = 5000; % mL (Peripheral compartment volume)
p.kc1 = .472; % hr-1 (clearance rate constant from central compartment)
p.kc2 = 0; % hr-1    (clearance rate constant from peripheral compartment)
p.k12 = .319;  % hr-1 (transport rate constant from central to peripheral)
p.k21 = .499; % hr-1  (transport rate constant from peripheral to central)
p.ka =  .764; % hr-1 (for tablet, use 1.27)  (absorption rate constant from gut to central)

% Initial conditions - softgel
D0 = 13.7*p.V1; % ug/ml * ml = ug (for tablet, use 12.9*p.V1)
% note - the initial dose amount in the paper was estimated as ug/ml and
%   proportional to central compartment volume; we're using this as an 
%   amount so we use CA0*V1
y0 = [0 0 0 D0]'; % y1-2 in ug/ml; y3-4 in ug

%% RUN SIMULATION - SOFTGEL
options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
[T1,Y1] = ode45(@m2C_Acet_eqns,[0:1/60:10],y0,options,p);

%% MASS BALANCE
% equations 1 and 2 produce concentrations, so multiply by the approporiate 
% volumes to get amounts
CurrentDrug(:,1) = Y1(:,1)*p.V1;
CurrentDrug(:,2) = Y1(:,2)*p.V2;
InitialDrug = y0(1)*p.V1 + y0(2)*p.V2 ;

% equations 3 and 4 are already in amount units, not concentration units; 
% so they don't need to be multiplied by volume
CurrentDrug(:,3) = Y1(:,3);
CurrentDrug(:,4) = Y1(:,4);

DrugIn = p.q*T1 + D0 ; % cumulative drug into system
DrugOut = CurrentDrug(:,3) ; % cumulative drug eliminated from system
BalanceD = DrugIn - DrugOut - CurrentDrug(:,1) - CurrentDrug(:,2) - CurrentDrug(:,4) + InitialDrug ; %(zero = balance)

%% CALCULATE AUC
% Trapezoidal rule
AUCoral1 = 0;
AUCoral2 = 0;
for i=1:(length(Y1)-1)
    AUCoral1 = AUCoral1 + 0.5*(Y1(i,1)+Y1(i+1,1))*(T1(i+1)-T1(i));
    AUCoral2 = AUCoral2 + 0.5*(Y1(i,2)+Y1(i+1,2))*(T1(i+1)-T1(i));
end

AUCoral1
AUCoral2
% no colon used here, so the results will output in the command window

%% VISUALIZE RESULTS - Softgel

% Figure 1

fig1 = figure;
ax1=subplot(2,2,1);
plot(ax1,T1,Y1(:,1),'k',T1,Y1(:,2),'r',T1,Y1(:,4)/p.V1,'b','linewidth',3)
title(ax1,'Drug Concentration')
ylabel(ax1,'[D] (ug/ml)')
xlabel(ax1,'time (hrs)')
lgd = legend('central', 'peripheral', 'gut');
lgd.Location = 'best';
lgd.Title.String = ['Compartment'];

ax2=subplot(2,2,2);
plot(ax2,T1,CurrentDrug(:,1),'k',T1,CurrentDrug(:,2),'r',T1,CurrentDrug(:,4),'b','linewidth',3)
title(ax2,'Total Drug Amount')
ylabel(ax2,'Total Drug (ug)')
xlabel(ax2,'time (hrs)')
lgd = legend('central', 'peripheral', 'gut');
lgd.Location = 'best';
lgd.Title.String = ['Compartment'];

ax3=subplot(2,2,3);
plot(ax3,T1,DrugIn,'g',T1,DrugOut,'c','linewidth',3)
title(ax3,'Cumulative Input/Output') 
ylabel(ax3,['Cumulative Drug in/out' newline 'of system (ug)'])
xlabel(ax3,'time (hrs)')
lgd = legend('input', 'output');
lgd.Location = 'best';
lgd.Title.String = ['Drug'];

ax4=subplot(2,2,4);
plot(ax4,T1,BalanceD,'m','linewidth',3)
title(ax4,'Molecular Balance') %(zero = balance)
ylabel(ax4,'Balance of Drug (ug)')
xlabel(ax4,'time (hrs)')





%% PARAMETERS - softgel-equivalent IV bolus

% No need to update the values of p - they will all be the same as for the
% softgel formulation; only the location of the bolus is different. Even
% the absorption rate ka can be left alone, because the concentration 
% in the gut is zero, so no absorption occurs anway.
%
% p.q = 0; % ug/hr
% p.V1 = 5000; % mL
% p.V2 = 5000; % mL
% p.kc1 = .472; % hr-1
% p.kc2 = 0; % hr-1
% p.k12 = .319;  % hr-1
% p.k21 = .499; % hr-1
% p.ka =  .764; % hr-1 

% Initial conditions - softgel-equivalent iv bolus
% Notice how the bolus is the same amount, but now is beginning as a
% concentration in the central compartment (1st element of y0) instead of 
% as an amount in the gut (4th element of y0).
%
% D0 = 13.7*p.V1; % ug/ml * ml = ug
y0 = [D0/p.V1 0 0 0]'; % y1-2 in ug/ml; y3-4 in ug

%% RUN SIMULATION - IV
options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
[T2,Y2] = ode45(@m2C_Acet_eqns,[0:1/60:10],y0,options,p);

%% MASS BALANCE
CurrentDrug(:,1) = Y2(:,1)*p.V1;
CurrentDrug(:,2) = Y2(:,2)*p.V2;
InitialDrug = y0(1)*p.V1 + y0(2)*p.V2 ;

CurrentDrug(:,3) = Y2(:,3);
CurrentDrug(:,4) = Y2(:,4);

DrugIn = p.q*T1  ; % cumulative drug into system
DrugOut = CurrentDrug(:,3) ; % cumulative drug eliminated from system
BalanceD = DrugIn - DrugOut - CurrentDrug(:,1) - CurrentDrug(:,2) - CurrentDrug(:,4) + InitialDrug ; %(zero = balance)

%% CALCULATE AUC
% Trapezoidal rule
AUCiv1 = 0;
AUCiv2 = 0;
for i=1:(length(Y2)-1)
    AUCiv1 = AUCiv1 + 0.5*(Y2(i,1)+Y2(i+1,1))*(T2(i+1,1)-T2(i,1));
    AUCiv2 = AUCiv2 + 0.5*(Y2(i,2)+Y2(i+1,2))*(T2(i+1,1)-T2(i,1));
end

AUCiv1
AUCiv2

%% CALCULATE BIOAVAILABILITY
Bioavailability = AUCoral1/AUCiv1
Bioavailability2 = AUCoral2/AUCiv2
% no colon used here, so the results will output in the command window

%% VISUALIZE RESULTS - IV

% Figure 2

fig2 = figure;
ax1=subplot(2,2,1);
plot(ax1,T2,Y2(:,1),'k',T2,Y2(:,2),'r',T2,Y2(:,4)/p.V1,'b','linewidth',3)
title(ax1,'Drug Concentration')
ylabel(ax1,'[D] (ug/ml)')
xlabel(ax1,'time (hrs)')
lgd = legend('central', 'peripheral', 'gut');
lgd.Location = 'best';
lgd.Title.String = ['Compartment'];

ax2=subplot(2,2,2);
plot(ax2,T2,CurrentDrug(:,1),'k',T2,CurrentDrug(:,2),'r',T2,CurrentDrug(:,4),'b','linewidth',3)
title(ax2,'Total Drug Amount')
ylabel(ax2,'Total Drug (ug)')
xlabel(ax2,'time (hrs)')
lgd = legend('central', 'peripheral', 'gut');
lgd.Location = 'best';
lgd.Title.String = ['Compartment'];

ax3=subplot(2,2,3);
plot(ax3,T2,DrugIn,'g',T2,DrugOut,'c','linewidth',3)
title(ax3,'Cumulative Input/Output') 
ylabel(ax3,['Cumulative Drug in/out' newline 'of system (ug)'])
xlabel(ax3,'time (hrs)')
lgd = legend('input', 'output');
lgd.Location = 'best';
lgd.Title.String = ['Drug'];

ax4=subplot(2,2,4);
plot(ax4,T2,BalanceD,'m','linewidth',3)
title(ax4,'Molecular Balance') %(zero = balance)
ylabel(ax4,'Balance of Drug (ug)')
xlabel(ax4,'time (hrs)')

%% VISUALIZE RESULTS - Compare routes of administration

% Figure 3

fig3 = figure;
ax9=subplot(1,2,1);
plot(ax9,T1,Y1(:,1),'k',T2,Y2(:,1),'r','linewidth',3)
title(ax9,'Drug Concentration (central)')
ylabel(ax9,'[D] (ug/ml)')
xlabel(ax9,'time (hrs)')
lgd = legend('Softgel', 'IV bolus');
lgd.Location = 'best';
lgd.Title.String = ['Formulation'];

ax10=subplot(1,2,2);
plot(ax10,T1,Y1(:,2),'k',T2,Y2(:,2),'r','linewidth',3)
title(ax10,'Drug Concentration (peripheral)')
ylabel(ax10,'[D] (ug/ml)')
xlabel(ax10,'time (hrs)')
lgd = legend('Softgel', 'IV bolus');
lgd.Location = 'best';
lgd.Title.String = ['Formulation'];


%% EXPORT VISUALIZATIONS
set(fig1,'Position',[0 0 600 450])
exportgraphics(fig1, "Fig1ba_Softgel.png",'Resolution',300); 

set(fig2,'Position',[0 0 600 450])
exportgraphics(fig2, "Fig2ba_IVbolus.png",'Resolution',300); 

set(fig3,'Position',[0 0 600 450])
exportgraphics(fig3, "Fig3ba_SoftgelIVComparison.png",'Resolution',300); 