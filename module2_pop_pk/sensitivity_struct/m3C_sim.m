function [outAUC,outT,outY] = m3C_sim(pin,OutputVar,TimeLen,VisFlag); 

% A large-molecule drug administered as a bolus dose D0 (fast injection)
% intravenously. The drug doesn't enter red blood cells, and undergoes
% first-order clearance (100 mL/hr) from the blood plasma (assume plasma is
% 2.5 L and makes up 50% of the blood). The drug can distribute from the
% bloodstream into a tumor of total size 0.5 L, 
% but 80% of the tumor is tumor cells or other cells, 
% and the remainder is extracellular space. Drug can also distribute from the
% bloodstream into the rest of the body (40 L, of which 88% consists of
% cells, and the remainder is extracellular space) with transport rate
% constants of 0.1 hr-1 (tumor to blood) and 0.5 hr-1 (rest of body to blood). 
% The reverse transport rate constants (i.e. from the bloodstream to the 
% peripheral compartments to the bloodstream have the same values adjusted 
% by a volume ratio (k12 = k21 * V2/V1)
% The drug is not cleared from the tissue
% compartments (i.e. the tumor and the rest of the body), only from the
% blood. In the blood, there is a protein that binds the drug with an
% equilibrium constant (Kd) of 6.5 nM and an unbinding (off) rate constant
% of 0.1 hr-1. The plasma protein is very large, and neither the plasma
% protein nor the bound form of the drug can leave the bloodstream (by
% clearance or by distribution). The plasma protein itself is initially
% present at 50 nM in the blood, and does not undergo clearance.

ylogflag = 0; % 1 for semilog y-axis

%% SIMULATION

y0 = [pin.D0/pin.V1 0 0 0 pin.B0 0]'; % nM - Initial conditions vector
% 1 = drug in blood
% 2 = drug in tumor 1
% 3 = drug in body
% 4 = drug in degr
% 5 = binding protein in blood
% 6 = drug-binding protein complex in blood

options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
[T,Y] = ode45(@m3C_eqns,[0:TimeLen/1000:TimeLen],y0,options,pin);

DrugInBolus = ones(length(T),1)*(pin.D0/pin.V1) ;

TotalFreeD(:,1) = Y(:,1)*pin.V1;
TotalFreeD(:,2) = Y(:,2)*pin.V2;
TotalFreeD(:,3) = Y(:,3)*pin.V3;
TotalBoundD(:,1)= Y(:,6)*pin.V1;

%% MASS BALANCE

DrugIn = pin.q*T + DrugInBolus*pin.V1 ; % cumulative drug into system
DrugOut = Y(:,4) ; % cumulative drug eliminated from system
BalanceD = DrugIn - DrugOut ...
                - TotalFreeD(:,1) - TotalFreeD(:,2) - TotalFreeD(:,3)  ...
                - TotalBoundD(:,1) ; %(zero = balance)


outAUC=trapz(T,Y(:,OutputVar));
outT=T;
outY=Y(:,OutputVar);


%% VISUALIZE RESULTS - 1

fig1=figure('visible',VisFlag);
ax1=subplot(4,2,1);
plot(ax1,T,Y(:,1),'k','linewidth',3)
title(ax1,'Concentration of Free Drug in Blood')
ylabel(ax1,'[D] (nM)')
xlabel(ax1,'time (hrs)')
if(ylogflag==1)
    set(gca,'yscale','log')
end

ax2=subplot(4,2,2);
plot(ax2,T,Y(:,1)*pin.V1+Y(:,6)*pin.V1,'k','linewidth',3)
title(ax2,'Total Amount of Drug in Blood')
ylabel(ax2,'Total Drug (nmol)')
xlabel(ax2,'time (hrs)')
if(ylogflag==1)
    set(gca,'yscale','log')
end

ax3=subplot(4,2,3);
plot(ax3,T,Y(:,2),'k','linewidth',3)
title(ax3,'Concentration of Free Drug in Tumor')
ylabel(ax3,'[D] (nM)')
xlabel(ax3,'time (hrs)')
if(ylogflag==1)
    set(gca,'yscale','log')
end

ax4=subplot(4,2,4);
plot(ax4,T,Y(:,2)*pin.V2,'k','linewidth',3)
title(ax4,'Total Amount of Drug in Tumor')
ylabel(ax4,'Total Drug (nmol)')
xlabel(ax4,'time (hrs)')
if(ylogflag==1)
    set(gca,'yscale','log')
end

ax7=subplot(4,2,5);
plot(ax7,T,Y(:,3),'k','linewidth',3)
title(ax7,'Concentration of Free Drug in Rest of Body')
ylabel(ax7,'[D] (nM)')
xlabel(ax7,'time (hrs)')
if(ylogflag==1)
    set(gca,'yscale','log')
end

ax8=subplot(4,2,6);
plot(ax8,T,Y(:,3)*pin.V3,'k','linewidth',3)
title(ax8,'Total Amount of Drug in Rest of Body')
ylabel(ax8,'Total Drug (nmol)')
xlabel(ax8,'time (hrs)')
if(ylogflag==1)
    set(gca,'yscale','log')
end

ax9=subplot(4,2,[7 8]);
plot(ax9,T,BalanceD,'k','linewidth',3)
title(ax9,'Molecular Balance') %(zero = balance)
ylabel(ax9,'Balance of Drug (nmol)')
xlabel(ax9,'time (hrs)')
if(ylogflag==1)
    set(gca,'yscale','log')
end

%% VISUALIZE RESULTS - 2 

fig2=figure('visible',VisFlag);

ax1=subplot(3,2,1);
plot(ax1,T,Y(:,1),'k','linewidth',3)
title(ax1,'Concentration of Free Drug in Blood')
ylabel(ax1,'[D] (nM)')
xlabel(ax1,'time (hrs)')
if(ylogflag==1)
    set(gca,'yscale','log')
end

ax2=subplot(3,2,2);
plot(ax2,T,Y(:,1)+Y(:,6),'k','linewidth',3)
title(ax2,'Concentration of Total Drug in Blood')
ylabel(ax2,'[D] + [DB] (nM)')
xlabel(ax2,'time (hrs)')
if(ylogflag==1)
    set(gca,'yscale','log')
end

ax3=subplot(3,2,3);
plot(ax3,T,Y(:,2),'k','linewidth',3)
title(ax3,'Concentration of Free Drug in Tumor')
ylabel(ax3,'[D] (nM)')
xlabel(ax3,'time (hrs)')
if(ylogflag==1)
    set(gca,'yscale','log')
end

ax7=subplot(3,2,4);
plot(ax7,T,Y(:,3),'k','linewidth',3)
title(ax7,'Concentration of Free Drug in Rest of Body')
ylabel(ax7,'[D] (nM)')
xlabel(ax7,'time (hrs)')
if(ylogflag==1)
    set(gca,'yscale','log')
end

ax9=subplot(3,2,[5 6]);
plot(ax9,T,BalanceD,'k','linewidth',3)
title(ax9,'Molecular Balance') %(zero = balance)
ylabel(ax9,'Balance of Drug (nmol)')
xlabel(ax9,'time (hrs)')
if(ylogflag==1)
    set(gca,'yscale','log')
end


%% VISUALIZE RESULTS - remaining stuff

fig3=figure('visible',VisFlag);
ax1=subplot(1,2,1);
plot(ax1,T,Y(:,6),'k','linewidth',3)
title(ax1,'Concentration of Boound Drug in Blood')
ylabel(ax1,'[DB] (nM)')
xlabel(ax1,'time (hrs)')
if(ylogflag==1)
    set(gca,'yscale','log')
end

ax2=subplot(1,2,2);
plot(ax2,T,Y(:,5),'k','linewidth',3)
title(ax2,'Concentration of Plasma Protein in Blood')
ylabel(ax2,'[B] (nM)')
xlabel(ax2,'time (hrs)')
if(ylogflag==1)
    set(gca,'yscale','log')
end

% output to mat files for R

%     save HW2Q2out.mat T Y BalanceD 
%     save HW2Q2out.mat T Y BalanceD

