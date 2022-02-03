clear all;

q = 1; % nmol/hr (continuous infusion)
D0 = 0; % nmol (bolus)
V = 1; % L
k = .1; % hr-1
y0 = [D0/V 0]'; % nM
p = [q V k]';

options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
[T1,Y1] = ode45(@m1C_eqns_withDegr,[0 10],y0,options,p);
DrugIn1 = q*T1 ; % cumulative drug into system
TotalD1 = Y1*V ; % Total drug in compartments
DrugOut1 = TotalD1(:,2) ; % cumulative drug eliminated from system
BalanceD1 = DrugIn1 - DrugOut1 - TotalD1(:,1) + y0(1)*V  ; %(zero = balance)
BalanceD1norm = BalanceD1./(DrugIn1+y0(1)*V) ; %(zero = balance)

q = 2;
V = 1; % L
D0 = 0; % nmol (bolus)
y0 = [D0/V 0]'; % nM
k = .1; % hr-1
p = [q V k]';
[T2,Y2] = ode45(@m1C_eqns_withDegr,[0 10],y0,options,p);
DrugIn2 = q*T2 ; % cumulative drug into system
TotalD2 = Y2*V ;
DrugOut2 = TotalD2(:,2) ; 
BalanceD2 = DrugIn2 - DrugOut2 - TotalD2(:,1) + y0(1)*V ; 
BalanceD2norm = BalanceD2./(DrugIn2+y0(1)*V) ; %(zero = balance)

q = 3;
V = 1; % L
D0 = 0; % nmol (bolus)
y0 = [D0/V 0]'; % nM
k = .1; % hr-1
p = [q V k]';
[T3,Y3] = ode45(@m1C_eqns_withDegr,[0 10],y0,options,p);
DrugIn3 = q*T3 ; % cumulative drug into system
TotalD3 = Y3*V ;
DrugOut3 = TotalD3(:,2) ; 
BalanceD3 = DrugIn3 - DrugOut3 - TotalD3(:,1) + y0(1)*V ; 
BalanceD3norm = BalanceD3./(DrugIn3+y0(1)*V) ; %(zero = balance)

figure;
ax1=subplot(3,2,1);
plot(ax1,T1,Y1(:,1),'k',T2,Y2(:,1),'b',T3,Y3(:,1),'r','linewidth',3)
title(ax1,'Concentration of Drug in Compartment (3 scenarios)')
ylabel(ax1,'[D] (nM)')
xlabel(ax1,'time (hrs)')

ax2=subplot(3,2,2);
plot(ax2,T1,TotalD1(:,1),'k',T2,TotalD2(:,1),'b',T3,TotalD3(:,1),'r','linewidth',3)
title(ax2,'Total Amount of Drug in System')
ylabel(ax2,'Total Drug (nmol)')
xlabel(ax2,'time (hrs)')

ax3=subplot(3,2,3);
plot(ax3,T1,DrugIn1,'k-',T2,DrugIn2,'b-',T3,DrugIn3,'r-',T1,DrugOut1,'k-.',T2,DrugOut2,'b-.',T3,DrugOut3,'r-.','linewidth',3)
title(ax3,'Cumulative Input (---)/Output(- -)') 
ylabel(ax3,'Cumulative Drug in/out of system (nmol)')
xlabel(ax3,'time (hrs)')

ax4=subplot(3,2,4);
plot(ax4,T1,BalanceD1,'k-',T2,BalanceD2,'r-',T3,BalanceD3,'b-','linewidth',3)
title(ax4,'Molecular Balance') %(zero = balance)
ylabel(ax4,'Balance of Drug (nmol)')
xlabel(ax4,'time (hrs)')

ax5=subplot(3,2,[5 6]);
plot(ax5,T1,BalanceD1norm,'k-',T2,BalanceD2norm,'r-',T3,BalanceD3norm,'b-','linewidth',3)
title(ax5,'Molecular Balance') %(zero = balance)
ylabel(ax5,'Balance of Drug (frac of drug in)')
xlabel(ax5,'time (hrs)')
