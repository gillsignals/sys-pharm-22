close all;
clear all;
% Driver file that calls an acetaminophen model simulation code for various
% specific scenarios:

runsim = 1; 
% Set this value to run different cases
% 1 = run two simulations (softgel, tablet) 
% 2 = run four simulations (softgel, tablet & iv equivalents for bioavailability) 
% 3 = run many values of kA - loop example
% 4 = run many values of CA0 - 2nd loop example

VisibleFigs = 1;
% Set this value to open each run's figure windows in the Acet_sim function
% (may not want to open figures if there are a lot!)
% 0 = figures not visible
% 1 = figures visible

SaveFigs = 1;
% Set this value to save each run's figure windows in the Acet_sim function 
% (may not want to save figures if there are a lot!)
% 0 = figures not visible
% 1 = figures visible

OralOrIV = 1; 
% This value determines delivery via oral or IV route
% 1 = oral
% 2 = IV

switch runsim
    case 1
        i=1; % Softgel formulation
        q0 = 0; % ug/hr
        V1x = 5000; % mL
        V2x = 5000; % mL
        CA0 = 13.7; % ug/ml (will be * ml = ug)
        kA =  .764; % hr-1
        
        OralOrIV = 1; % 1 = oral; 2 = IV
        [AUC1(i),AUC2(i)] = Acet_sim(q0,CA0,kA,V1x,V2x,OralOrIV,VisibleFigs,SaveFigs,"Softgel");

        i=2; % Tablet formulation
%       These values are unchanged for tablet formulation    
%         q0 = 0; % ug/hr
%         V1x = 5000; % mL
%         V2x = 5000; % mi
        CA0 = 12.9; % ug/ml * ml = ug
        kA =  1.27; % hr-1

        OralOrIV = 1; % 1 = oral; 2 = IV
        [AUC1(i),AUC2(i)] = Acet_sim(q0,CA0,kA,V1x,V2x,OralOrIV,VisibleFigs,SaveFigs,"Tablet");
        
        AUC1
        AUC2
        
        fig2 = figure;
        subplot(1,2,1);
        bar(AUC1);
        title(gca,'Central Compartment AUC') 
        ylabel(gca,'AUC (ug*hr/ml)')
        xlabel(gca,'Scenario')
        set(gca,'XTickLabel',{'Softgel','Tablet'});

        subplot(1,2,2);
        bar(AUC2);
        title(gca,'Peripheral Compartment AUC') 
        ylabel(gca,'AUC (ug*hr/ml)')
        xlabel(gca,'Scenario')
        set(gca,'XTickLabel',{'Softgel','Tablet'});
        
        % EXPORT VISUALIZATION
        set(fig2,'Position',[0 0 600 450])
        exportgraphics(fig2, 'Fig2_CentralPeripheralAUC.png','Resolution',300); 

        
    case 2
        i=1;
        q0 = 0; % ug/hr
        V1x = 5000; % mL
        V2x = 5000; % mL
        CA0 = 13.7; % ug/ml (will be * ml = ug)
        kA =  .764; % hr-1

        OralOrIV = 1; % 1 = oral; 2 = IV
        [AUC1(i),AUC2(i)] = Acet_sim(q0,CA0,kA,V1x,V2x,OralOrIV,VisibleFigs,SaveFigs,"Softgel_Oral");

        i=2;        
        OralOrIV = 2; % 1 = oral; 2 = IV
        [AUC1(i),AUC2(i)] = Acet_sim(q0,CA0,kA,V1x,V2x,OralOrIV,VisibleFigs,SaveFigs,"SoftgelEquiv_IV");

        i=3;
%         q0 = 0; % ug/hr
%         V1x = 5000; % mL
%         V2x = 5000; % mL
        CA0 = 12.9; % ug/ml * ml = ug
        kA =  1.27; % hr-1

        OralOrIV = 1; % 1 = oral; 2 = IV
        [AUC1(i),AUC2(i)] = Acet_sim(q0,CA0,kA,V1x,V2x,OralOrIV,VisibleFigs,SaveFigs,"Tablet_Oral");
                
        i=4;
        OralOrIV = 2; % 1 = oral; 2 = IV
        [AUC1(i),AUC2(i)] = Acet_sim(q0,CA0,kA,V1x,V2x,OralOrIV,VisibleFigs,SaveFigs,"TabletEquiv_IV");
        
        AUC1
        AUC2
        
        fig2 = figure;
        subplot(1,2,1);
        bar(AUC1);
        title(gca,'Central Compartment AUC') 
        ylabel(gca,'AUC (ug*hr/ml)')
        xlabel(gca,'Scenario')
        set(gca,'XTickLabel',{'Softgel','Softgel equiv iv','Tablet','Tablet equiv iv'});

        subplot(1,2,2);
        bar(AUC2);
        title(gca,'Peripheral Compartment AUC') 
        ylabel(gca,'AUC (ug*hr/ml)')
        xlabel(gca,'Scenario')
        set(gca,'XTickLabel',{'Softgel','Softgel equiv iv','Tablet','Tablet equiv iv'});

        bioavail_softgel_central = AUC1(1)/AUC1(2)
        bioavail_tablet_central = AUC1(3)/AUC1(4)
        bioavail_softgel_peripheral = AUC2(1)/AUC2(2)
        bioavail_tablet_peripheral = AUC2(3)/AUC2(4)

        % EXPORT VISUALIZATION
        set(fig2,'Position',[0 0 600 450])
        exportgraphics(fig2, 'Fig2_CentralPeripheralAUC.png','Resolution',300); 

    case 3

        q0 = 0; % ug/hr
        V1x = 5000; % mL
        V2x = 5000; % mL
        CA0 = 12.9; % ug/ml * ml = ug
        kA = linspace(0.5,2.5,5);

        VisibleFigs = 0; % set to 0 if we don't want to open figure windows every loop   
%         SaveFigs = 0;    % set to 0 if we don't want to save figures every loop   

        for i=1:length(kA)
            [AUC1(i),AUC2(i)] = Acet_sim(q0,CA0,kA(i),V1x,V2x,OralOrIV,VisibleFigs,SaveFigs,strcat('kA_run',num2str(i)));
        end
        
        AUC1
        AUC2
        
        fig2 = figure;
        plot(kA,AUC1,'k');
        hold on;
        scatter(kA,AUC1,'r');
        title(gca,'AUC vs kA') 
        ylabel(gca,'AUC (ug*hr/ml)')
        xlabel(gca,'kA (hr^-^1)')

        % EXPORT VISUALIZATION
        set(fig2,'Position',[0 0 600 450])
        exportgraphics(fig2, 'Fig2_AUCvskA.png','Resolution',300); 

	case 4

        q0 = 0; % ug/hr
        V1x = 5000; % mL
        V2x = 5000; % mL
        CA0 = logspace(0,2,5); % ug/ml * ml = ug
        kA = .764;

        VisibleFigs = 0; % set to 0 if we don't want to open figure windows every loop   
        SaveFigs = 0;    % set to 0 if we don't want to save figures every loop   

        for i=1:length(CA0)
            [AUC1(i),AUC2(i)] = Acet_sim(q0,CA0(i),kA,V1x,V2x,OralOrIV,VisibleFigs,SaveFigs,strcat('CA0_run',num2str(i)));
        end
        
        AUC1
        AUC2
        
        fig2 = figure;
        ax1=subplot(1,2,1);
        plot(ax1,CA0,AUC1,'k');
        hold on;
        scatter(ax1,CA0,AUC1,'r');
        title(ax1,'AUC vs Dose; linear x-axis') 
        ylabel(ax1,'AUC (ug*hr/ml)')
        xlabel(ax1,'Drug Dose (ug/ml)')

        ax2=subplot(1,2,2);
        plot(ax2,CA0,AUC1,'k');
        hold on;
        scatter(ax2,CA0,AUC1,'r');
        set(ax2,'XScale','log')
        title(ax2,'AUC vs Dose; log x-axis') 
        ylabel(ax2,'AUC (ug*hr/ml)')
        xlabel(ax2,'Drug Dose (ug/ml)')
        
        % EXPORT VISUALIZATION
        set(fig2,'Position',[0 0 600 450])
        exportgraphics(fig2, 'Fig2_AUCvsDose.png','Resolution',300); 

end
