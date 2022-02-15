close all;
clear all;

%% DEFINE CASES/FLAGS

CaseFlag = 1;
  % CASE ONE -- adult, cup of coffee
  % CASE TWO -- adult, cup of coffee (loop formulation)
  % CASE THREE -- plotting concentrations over time
  % CASE FOUR -- plotting concentrations over time (loop formulation)
  % CASE FIVE -- show real time
  % CASE SIX -- 2d plots - use one output
  % CASE ELEVEN-THIRTEEN - vary halflives (for ppt slides)
% outputs: 
    % cases 1,2,6: Max Concentration
    % cases 3,4: concentration profile
    % case  5: timesteps & concentration profile

VisFlag = 'off'; % off = suppress simulation figures; on = don't suppress

TimeLen = 2; % hrs

%% RUN CASES

switch CaseFlag 

case 1
%% CASE ONE -- adult, cup of coffee
%========================
% OUTPUT: peak concentration of caffeine
% VARY: clearance halflife

% Parameters
  weight = 70 ; %kg
  Vd = 0.6 * weight ; % L

  abshalflife = 7/60; % 7 mins in hr (ABSORPTION HALFLIFE)
  halflife = 5; % hr (CLEARANCE HALFLIFE)
  dose = 310; % mg - bolus dose, grande coffee

% Run simulations
  [a,~] = Sim_Caffeine(abshalflife,Vd,halflife*0.01,dose,CaseFlag,VisFlag,TimeLen);
  [b,~] = Sim_Caffeine(abshalflife,Vd,halflife*0.1,dose,CaseFlag,VisFlag,TimeLen);
  [c,~] = Sim_Caffeine(abshalflife,Vd,halflife*1,dose,CaseFlag,VisFlag,TimeLen);
  [d,~] = Sim_Caffeine(abshalflife,Vd,halflife*10,dose,CaseFlag,VisFlag,TimeLen);
  [e,~] = Sim_Caffeine(abshalflife,Vd,halflife*100,dose,CaseFlag,VisFlag,TimeLen);

% Plot results from all simulations  
  toplot = [a b c d e];
  figure;
  bar(toplot);
  ax = gca; % assign a handle (ax) for the current axes
  ax.FontSize = 14;
  xticklabels({halflife*0.01,halflife*0.1,halflife,halflife*10,halflife*100});
  ylabel('Peak caffeine in bloodstream (mg/L)','FontSize',16)
  xlabel('Clearance Halflife (hrs)','FontSize',16)

%========================


case 2
%% CASE TWO -- adult, cup of coffee (loop formulation)

%========================
% OUTPUT: peak concentration of caffeine
% VARY: clearance halflife

% Parameters
weight = 70 ; %kg
Vd = 0.6 * weight ; % L

abshalflife = 7/60; % 7 mins in hr (ABSORPTION HALFLIFE)
halflife = 5; % hr (CLEARANCE HALFLIFE)
dose = 310; % mg - bolus dose, grande coffee

% Different ways of defining the input parameters
% a = [0.01,0.1,1,10,100]*halflife ;
% a = logspace(-2,2,50)*halflife;
a = [2,3,4,5,6,7] ;
% a = linspace(1,10,50);

% Run simulations
for i=1:length(a)
    [b(i),~] = Sim_Caffeine(abshalflife,Vd,a(i),dose,CaseFlag,VisFlag,TimeLen);
end

% Plot results from all simulations
figure;
bar(b);
ax = gca; % assign a handle (ax) for the current axes
ax.FontSize = 14;
xticklabels(a);
ylabel('Peak caffeine in bloodstream (mg/L)','FontSize',16)
xlabel('Clearance halflife (hr)','FontSize',16)

figure;
scatter(a,b,'k');
set(gca,'xscale','log')
ax = gca; % assign a handle (ax) for the current axes
ax.FontSize = 14;
ylabel('Peak caffeine in bloodstream (mg/L)','FontSize',16)
xlabel('Clearance halflife (hr)','FontSize',16)

%========================


case 3
%% CASE THREE -- plotting concentrations over time
%========================
% OUTPUT: concentration profiles of caffeine over time
% VARY: clearance halflife

% Parameters
weight = 70 ; %kg
Vd = 0.6 * weight ; % L

abshalflife = 7/60; % 7 mins in hr (ABSORPTION HALFLIFE)
halflife = 5; % hr (CLEARANCE HALFLIFE)
dose = 310; % mg - bolus dose, grande coffee

% Run simulations
[a,~] = Sim_Caffeine(abshalflife,Vd,halflife*0.01,dose,CaseFlag,VisFlag,TimeLen);
[b,~] = Sim_Caffeine(abshalflife,Vd,halflife*0.1,dose,CaseFlag,VisFlag,TimeLen);
[c,~] = Sim_Caffeine(abshalflife,Vd,halflife*1,dose,CaseFlag,VisFlag,TimeLen);
[d,~] = Sim_Caffeine(abshalflife,Vd,halflife*10,dose,CaseFlag,VisFlag,TimeLen);
[e,~] = Sim_Caffeine(abshalflife,Vd,halflife*100,dose,CaseFlag,VisFlag,TimeLen);

% Plot results from all simulations  
figure;
plot(a,'LineWidth',3);
hold on;
plot(b,'LineWidth',3);
plot(c,'LineWidth',3);
plot(d,'LineWidth',3);
plot(e,'LineWidth',3);
ax = gca; % assign a handle (ax) for the current axes
ax.FontSize = 14;
ylabel('[Caffeine in bloodstream] (mg/L)','FontSize',16)
xlabel('Time (timesteps)','FontSize',16)
lgd = legend('0.01','0.1','1','10','100');
lgd.Location = 'best';
lgd.Title.String = ['Clearance' newline 'halflife' newline '(multiple of' newline 'baseline)'];

% %========================

case 4
%% CASE FOUR -- plotting concentrations over time (loop formulation)
%========================
% OUTPUT: concentration profiles of caffeine over time
% VARY: clearance halflife

% Parameters
weight = 70 ; %kg
Vd = 0.6 * weight ; % L

abshalflife = 7/60; % 7 mins in hr (ABSORPTION HALFLIFE)
halflife = 5; % hr (CLEARANCE HALFLIFE)
dose = 310; % mg - bolus dose, grande coffee

% Parameters to scan: pick one 
% a = [0.01,0.1,1,10,20]*halflife;
% a = [2,3,4,5,6,7] ;
% a = logspace(-2,2,50)*halflife;
a = linspace(2,7,5);

% Run simulations
for i=1:length(a)
    [b(:,i),~] = Sim_Caffeine(abshalflife,Vd,a(i),dose,CaseFlag,VisFlag,TimeLen);
end

% Plot results from all simulations  
figure;
plot(b,'LineWidth',3);
ax = gca; % assign a handle (ax) for the current axes
ax.FontSize = 14;
ylabel('[Caffeine in bloodstream] (mg/L)','FontSize',16)
xlabel('Time (timesteps)','FontSize',16)
lgd = legend(string(a));
lgd.Location = 'best';
lgd.Title.String = ['Clearance' newline 'halflife' newline 't1/2 (hr)'];


%========================


case 5
%% CASE FIVE -- show real time
%========================
% OUTPUT: time, concentration profiles of caffeine over time
% VARY: clearance halflife

% Parameters
weight = 70 ; %kg
Vd = 0.6 * weight ; % L

abshalflife = 7/60; % 7 mins in hr (ABSORPTION HALFLIFE)
halflife = 5; % hr (CLEARANCE HALFLIFE)
dose = 310; % mg - bolus dose, grande coffee

% Run simulations
[t1,a] = Sim_Caffeine(abshalflife,Vd,halflife*0.01,dose,CaseFlag,VisFlag,TimeLen);
[t2,b] = Sim_Caffeine(abshalflife,Vd,halflife*0.1,dose,CaseFlag,VisFlag,TimeLen);
[t3,c] = Sim_Caffeine(abshalflife,Vd,halflife*1,dose,CaseFlag,VisFlag,TimeLen);
[t4,d] = Sim_Caffeine(abshalflife,Vd,halflife*10,dose,CaseFlag,VisFlag,TimeLen);
[t5,e] = Sim_Caffeine(abshalflife,Vd,halflife*100,dose,CaseFlag,VisFlag,TimeLen);

% Plot results from all simulations  
figure;
plot(t1,a, 'k','LineWidth',3);
hold on;
plot(t2,b,'LineWidth',3);
plot(t3,c,'LineWidth',3);
plot(t4,d,'LineWidth',3);
plot(t5,e, 'r','LineWidth',3);
ax = gca; % assign a handle (ax) for the current axes
ax.FontSize = 14;
ylabel('[Caffeine in bloodstream] (mg/L)','FontSize',16)
xlabel('Time (hrs)','FontSize',16)
lgd = legend('0.01','0.1','1','10','100');
lgd.Location = 'best';
lgd.Title.String = ['Clearance' newline 'halflife' newline '(multiple of' newline 'baseline)'];


% ========================

case 6
%% CASE SIX - 2d plots - use one output
%========================
% OUTPUT: peak concentration of caffeine
% VARY: volume, clearance halflife

% Parameters
weight = 70 ; %kg
Vd = 0.6 * weight ; % L

abshalflife = 7/60; % 7 mins in hr (ABSORPTION HALFLIFE)
halflife = 5; % hr (CLEARANCE HALFLIFE)
dose = 310; % mg - bolus dose, grande coffee

p = logspace(-1,1,5); % global sensitivity - two orders of magnitude

% Run simulations
for i=1:length(p)
for j=1:length(p)
    [a(i,j),~] = Sim_Caffeine(abshalflife,Vd*p(i),halflife*p(j),dose,CaseFlag,VisFlag,TimeLen);
end
end

% Plot results from all simulations  
figc = figure;
[M,c] = contourf(p*Vd,p*halflife,a');
colormap(parula)
colorbar;
set(gca,'xscale','log')
set(gca,'yscale','log')
ax = gca; % assign a handle (ax) for the current axes
ax.FontSize = 14;
clabel(M,c,'FontSize',14)
ylabel('clearance halflife (hrs)','FontSize',16)
xlabel('Volume (L)','FontSize',16)
title('Peak Caffeine Concentration (mg/L) - global sensitivity','FontSize',16)

saveas (figc,'2d_sensitivity_contour.png')

figc2 = figure;
h1 = heatmap(round(p*Vd,1),round(p*halflife,1),a');
h1.FontSize = 14
h1.Title = 'Peak Caffeine Concentration (mg/L) - global sensitivity';
h1.XLabel = 'Volume (L)';
h1.YLabel = 'clearance halflife (hrs)';

case 11
%% CASE ELEVEN -- show real time, vary absorption
%========================
% OUTPUT: time, concentration profiles of caffeine in gut over time
% VARY: clearance halflife

% Parameters
weight = 70 ; %kg
Vd = 0.6 * weight ; % L

abshalflife = 7/60; % 7 mins in hr (ABSORPTION HALFLIFE)
halflife = 5; % hr (CLEARANCE HALFLIFE)
dose = 310; % mg - bolus dose, grande coffee

% Run simulations
[t1,a] = Sim_Caffeine(abshalflife,Vd,halflife,dose,CaseFlag,VisFlag,TimeLen);
[t2,b] = Sim_Caffeine(abshalflife*2,Vd,halflife,dose,CaseFlag,VisFlag,TimeLen);
[t3,c] = Sim_Caffeine(abshalflife*3,Vd,halflife,dose,CaseFlag,VisFlag,TimeLen);

% Plot results from all simulations  
fig1 = figure; % give the figure a handle (fig1) to refer to it later
plot(t1*60,a, 'k','LineWidth',4); % note *60 to convert hrs to mins
hold on;
plot(t2*60,b, 'b','LineWidth',4);
plot(t3*60,c, 'r','LineWidth',4);
ax = gca; % assign a handle (ax) for the current axes
ax.FontSize = 14;
title('Amount of Caffeine in the gut','FontSize',16)
ylabel('Caffeine in gut (mg)','FontSize',16)
xlabel('time (mins)','FontSize',16)

saveas (fig1,'caffeine_gut_vary_absorption.png')

% ========================

case 12
%% CASE TWELVE -- show real time, vary absorption
%========================
% OUTPUT: time, concentration profiles of caffeine in gut over time
% VARY: clearance halflife

% Parameters
weight = 70 ; %kg
Vd = 0.6 * weight ; % L

abshalflife = 7/60; % 7 mins in hr (ABSORPTION HALFLIFE)
halflife = 5; % hr (CLEARANCE HALFLIFE)
dose = 310; % mg - bolus dose, grande coffee

% Run simulations
[t1,a] = Sim_Caffeine(abshalflife,Vd,halflife,dose,CaseFlag,VisFlag,TimeLen);
[t2,b] = Sim_Caffeine(abshalflife*2,Vd,halflife,dose,CaseFlag,VisFlag,TimeLen);
[t3,c] = Sim_Caffeine(abshalflife*3,Vd,halflife,dose,CaseFlag,VisFlag,TimeLen);

% Plot results from all simulations  
fig1 = figure; % give the figure a handle (fig1) to refer to it later
plot(t1*60,a, 'k','LineWidth',4); % note *60 to convert hrs to mins
hold on;
plot(t2*60,b, 'b','LineWidth',4);
plot(t3*60,c, 'r','LineWidth',4);
ax = gca; % assign a handle (ax) for the current axes
ax.FontSize = 14;
title('Concentration of Caffeine in the bloodstream','FontSize',16)
ylabel('Caffeine in bloodstream (mg/L)','FontSize',16)
xlabel('time (mins)','FontSize',16)

saveas (fig1,'caffeine_bloodstream_vary_absorption.png')

% ========================

case 13
%% CASE THIRTEEN -- show real time, vary clearance
%========================
% OUTPUT: time, concentration profiles of caffeine in gut over time
% VARY: clearance halflife

% Parameters
weight = 70 ; %kg
Vd = 0.6 * weight ; % L

abshalflife = 7/60; % 7 mins in hr (ABSORPTION HALFLIFE)
halflife = 5; % hr (CLEARANCE HALFLIFE)
dose = 310; % mg - bolus dose, grande coffee

% Run simulations
[t1,a] = Sim_Caffeine(abshalflife,Vd,halflife,dose,CaseFlag,VisFlag,TimeLen);
[t2,b] = Sim_Caffeine(abshalflife,Vd,halflife*2,dose,CaseFlag,VisFlag,TimeLen);
[t3,c] = Sim_Caffeine(abshalflife,Vd,halflife/2,dose,CaseFlag,VisFlag,TimeLen);

% Plot results from all simulations  
fig1 = figure; % give the figure a handle (fig1) to refer to it later
plot(t1*60,a, 'k','LineWidth',4); % note *60 to convert hrs to mins
hold on;
plot(t2*60,b, 'b','LineWidth',4);
plot(t3*60,c, 'r','LineWidth',4);
ax = gca; % assign a handle (ax) for the current axes
ax.FontSize = 14;
title('Concentration of Caffeine in the bloodstream','FontSize',16)
ylabel('Caffeine in bloodstream (mg/L)','FontSize',16)
xlabel('time (mins)','FontSize',16)

saveas (fig1,'caffeine_bloodstream_vary_clearance.png')

% ========================


end