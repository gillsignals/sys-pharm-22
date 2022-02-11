function [out1,out2] = Acet_sim(q0,CA0,kA,V1x,V2x,OralOrIV,OnOrOff,SaveFigs,SimName); 
% Simulate one 'run' (simulation) of the Acetaminophen model. 
%
% This function takes in key parameters, defines others, and calls the
% solver for the simulation.
%
% The function returns two outputs (the AUC for compartment 1 and 2); 
%  and receives five parameters, two flags, and a name string as inputs. 
%
% The line above can be modified as needed to build a function that 
% returns/takes in more, fewer, or different outputs or parameters. 
% As a general rule, pass the parameters/outputs you think you will need,
% rather than passing everything.
%
%% THIS FUNCTION RETURNS:
% out1 = AUC in central compartment
% out2 = AUC in peripheral compartment
%
%% ARGUMENTS
% q0 = infusion rate
% CA0 = initial concentration in gut compartment (or equivalent for IV)
% kA = absorption rate constant
% V1x = Volume of central compartment
% V2x = Volume of peripheral compartment
% OralOrIV = flag (0 for Oral delivery; 1 for IV delivery)
% OnOrOff = figure visibility (0 for invisible; 1 for visible)
% SaveFigs = exporting figures to files (0 for no; 1 for yes)
% SimName = short text string to identify this simulation; used in naming
%  output figure files

%% PARAMETER VALUES
% note some of these values are passed in from the driver
% it's not necessary to change the names (e.g. q = q0), you could
% just use the names above. We change the names only to emphasize
% where they come from, and to leave it flexible and easy to change which
% parameters get passed.

p.q = q0; % ug/hr
p.V1 = V1x; % mL
p.V2 = V2x; % mL
p.kc1 = .472; % hr-1
p.kc2 = 0; % hr-1
p.k12 = .319;  % hr-1
p.k21 = .499; % hr-1
p.ka =  kA; % hr-1

D0 = CA0*p.V1; % ug/ml * ml = ug

% Initial conditions
if OralOrIV == 1 % oral
    y0 = [0 0 0 D0]'; % y1,y2 in ug/ml; y3,y4 in ug
elseif OralOrIV == 2 % IV
    y0 = [D0/p.V1 0 0 0]'; % y1,y2 in ug/ml; y3,y4 in ug
end    

SimLength = 10; % hours
ReportingTimeStep = 1.0/60; % ask solver to return result every minute

%% CALL SOLVER

options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
[T1,Y1] = ode45(@Acet_eqns,[0:ReportingTimeStep:SimLength],y0,options,p);

%% USE RESULTS TO CALCULATE MASS BALANCE and AUC

% equations 1 and 2 produce concentrations, so multiply by the approporiate 
% volumes to get amounts
CurrentDrug(:,1) = Y1(:,1)*p.V1;
CurrentDrug(:,2) = Y1(:,2)*p.V2;
InitialDrug = y0(1)*p.V1 + y0(2)*p.V2 ;

% equations 3 and 4 are already in amount units, not concentration units; 
% so they don't need to be multiplied by volume
CurrentDrug(:,3) = Y1(:,3);
CurrentDrug(:,4) = Y1(:,4);

if OralOrIV == 1 % oral
    DrugIn = p.q*T1 + D0 ; % cumulative drug into system (continuous infusion + bolus)
elseif OralOrIV == 2 % IV
    DrugIn = p.q*T1  ; % D0 already accounted for in initial drug
end    
DrugOut = CurrentDrug(:,3) ; % cumulative drug eliminated from system
BalanceD = DrugIn - DrugOut - CurrentDrug(:,1) - CurrentDrug(:,2) - CurrentDrug(:,4) + InitialDrug ; %(zero = balance)

% calculate AUC by integrating the concentration curve (trapezoidal rule)
AUC1 = 0;
AUC2 = 0;
for i=1:(length(Y1)-1)
    AUC1 = AUC1 + 0.5*(Y1(i,1)+Y1(i+1,1))*(T1(i+1)-T1(i));
    AUC2 = AUC2 + 0.5*(Y1(i,2)+Y1(i+1,2))*(T1(i+1)-T1(i));
end

%% VISUALIZATION

if SaveFigs == 1 || OnOrOff == 1 % if both zero, don't bother making the figures at all!

if OnOrOff == 1
    fig1 = figure;
else
    fig1 = figure('visible','off');
end

% figure;

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

%% EXPORT VISUALIZATION
if SaveFigs == 1
    set(fig1,'Position',[0 0 600 450]);
    exportgraphics(fig1, strcat('Fig1_',SimName,'.png'),'Resolution',300); 
end

end

%% RETURN OUTPUTS

out1 = AUC1;
out2 = AUC2;

end
