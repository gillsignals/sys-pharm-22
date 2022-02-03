clear all;

q = 1; % nmol/hr
V = 1; % L
kc1 = .1; % hr-1
kc2 = kc1; % hr-1 - same clearance as A
% kc2 = kc1/5; % hr-1 - slower clearance for B
kab = 1; % M-1 hr-1
y0 = [0 0]'; % nM
p = [q V kc1 kc2 kab]';

options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
[T1,Y1] = ode45(@m1C_eqns_Metab,[0 10],y0,options,p);
TotalD1 = Y1*V ;

% q = 2;
% kc2 = .01;
kab = .1; % M-1 hr-1 - lower conversion rate
p = [q V kc1 kc2 kab]';
[T2,Y2] = ode45(@m1C_eqns_Metab,[0 10],y0,options,p);
TotalD2 = Y2*V ;

% q = 3;
% kc2 = .001;
kab = .01; % M-1 hr-1 - lower conversion rate
p = [q V kc1 kc2 kab]';
[T3,Y3] = ode45(@m1C_eqns_Metab,[0 10],y0,options,p);
TotalD3 = Y3*V ;

% figure;
% plot(T1,Y1(:,1),'k',T1,Y1(:,2),'b',T1,Y1(:,3),'r','linewidth',3)

%plot(T1,TotalD1(:,1),'k',T2,TotalD2(:,1),'b',T3,TotalD3(:,1),'r','linewidth',3)

figure;
ax1=subplot(2,2,1);
plot(ax1,T1,Y1(:,1),'k',T2,Y2(:,1),'b',T3,Y3(:,1),'r','linewidth',3)
title(ax1,'Concentration of A')
ylabel(ax1,'[D] (nM)')
xlabel(ax1,'time (hrs)')

ax2=subplot(2,2,2);
plot(ax2,T1,Y1(:,2),'k',T2,Y2(:,2),'b',T3,Y3(:,2),'r','linewidth',3)
title(ax2,'Concentration of B')
ylabel(ax2,'[D] (nM)')
xlabel(ax2,'time (hrs)')

ax3=subplot(2,2,3);
plot(ax3,T1,TotalD1(:,1),'k',T2,TotalD2(:,1),'b',T3,TotalD3(:,1),'r','linewidth',3)
title(ax3,'Total Amount of Free Drug (A) in Compartment')
ylabel(ax3,'Total Drug (nmol)')
xlabel(ax3,'time (hrs)')

ax4=subplot(2,2,4);
plot(ax4,T1,TotalD1(:,1)+TotalD1(:,2),'k',T2,TotalD2(:,1)+TotalD2(:,2),'b',T3,TotalD3(:,1)+TotalD3(:,2),'r','linewidth',3)
title(ax4,'Total Amount of Total Drug (A+B) in Compartment')
ylabel(ax4,'Total Drug (nmol)')
xlabel(ax4,'time (hrs)')