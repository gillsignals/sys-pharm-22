clear all;

%% PARAMETERS

OutputFlag = 1; 
% outputs: 
    % cases 1,2,6: Max Concentration
    % cases 3,4: concentration profile
    % case  5: timesteps & concentration profile
VisFlag = 'off';

TimeLen = 2; % hrs

%% DISTRIBUTIONS

numPatients = 100;

% uniform weight distribution
lb = 50;
ub = 90;

% normal weight distribution
mean = 70;
stddev= 20;

%% PARAMETERS

abshalflife = 7/60; % 7 mins in hr
halflife = 5; % hr
dose = 310; % mg - bolus dose, grande coffee


%% RUN SIMULATIONS

for k=1:numPatients
    weightu(k) = random('Uniform',lb,ub) ; % uniform
    weightn(k) = normrnd(mean,stddev) ; % normal
%     Vda(k) = 0.6*weightu(k) ;
    Vda(k) = 0.6*weightu(k) * normrnd(1,.2); % added interindividual variability
    a(k) = Sim_Caffeine(abshalflife,Vda(k),halflife,dose,OutputFlag,VisFlag,TimeLen);
%     Vdb(k) = 0.6*weightn(k) ;
    Vdb(k) = 0.6*weightn(k) * normrnd(1,.2); % added interindividual variability
    b(k) = Sim_Caffeine(abshalflife,Vdb(k),halflife,dose,OutputFlag,VisFlag,TimeLen);
end

%% VISUALIZATIONS -- Peak Concentrations

fig1 = figure;
scatter(weightn,b,'r') ;
hold on;
scatter(weightu,a,'k') ;
ax = gca; % assign a handle (ax) for the current axes
ax.FontSize = 14;
title(ax,'Peak Concentrations vs Weight (red = normal) (black = uniform)','FontSize',16) 
ylabel(ax,'Peak Caffeine Concentration (mg/L)','FontSize',16)
xlabel(ax,'Weight (kg)','FontSize',16)

fig2 = figure;
scatter(weightu,a,'k') ;
ax = gca; % assign a handle (ax) for the current axes
ax.FontSize = 14;
title(gca,'Peak Concentrations vs Weight (uniform)','FontSize',16) 
ylabel(gca,'Peak Caffeine Concentration (mg/L)','FontSize',16)
xlabel(gca,'Weight (kg)','FontSize',16)

fig3 = figure;
scatter(weightn,b,'r') ;
ax = gca; % assign a handle (ax) for the current axes
ax.FontSize = 14;
title(gca,'Peak Concentrations vs Weight (normal)','FontSize',16) 
ylabel(gca,'Peak Caffeine Concentration (mg/L)','FontSize',16)
xlabel(gca,'Weight (kg)','FontSize',16)

saveas (fig1,'PeakConcs_both.png')
saveas (fig2,'PeakConcs_uniform.png')
saveas (fig3,'PeakConcs_normal.png')

%% VISUALIZATIONS -- Volumes of distribution

fig4 = figure;
scatter(weightn,Vdb,'r') ;
hold on;
scatter(weightu,Vda,'k') ;
ax = gca; % assign a handle (ax) for the current axes
ax.FontSize = 14;
title(ax,'Volume of distribution vs Weight (red = normal) (black = uniform)','FontSize',16) 
ylabel(ax,'Volume of Distribution (L)','FontSize',16)
xlabel(ax,'Weight (kg)','FontSize',16)

fig5 = figure;
scatter(weightu,Vda,'k') ;
ax = gca; % assign a handle (ax) for the current axes
ax.FontSize = 14;
title(gca,'Volume of distribution vs Weight (uniform)','FontSize',16) 
ylabel(gca,'Volume of Distribution (L)','FontSize',16)
xlabel(gca,'Weight (kg)','FontSize',16)

fig6 = figure;
scatter(weightn,Vdb,'r') ;
ax = gca; % assign a handle (ax) for the current axes
ax.FontSize = 14;
title(gca,'Volume of distribution vs Weight (normal)','FontSize',16) 
ylabel(gca,'Volume of Distribution (L)','FontSize',16)
xlabel(gca,'Weight (kg)','FontSize',16)

saveas (fig4,'VolDist_both.png')
saveas (fig5,'VolDist_uniform.png')
saveas (fig6,'VolDist_normal.png')


%========================




