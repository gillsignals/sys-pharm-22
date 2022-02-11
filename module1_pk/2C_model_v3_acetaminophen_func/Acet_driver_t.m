clear all;

flag = 3; 
% set this value to run different cases
% 1 = run two simulations (softgel, tablet) 
% 2 = run two simulations (oral vs iv, softgel)
% 3 = run two simulations (oral vs iv, tablet)

VisibleFigs = 0;
% set this value to open figure windows in the Acet_sim function
% 0 = figures not visible
% 1 = figures visible

OralOrIV = 1; 
% this value determines delivery via oral or IV
% 1 = oral
% 2 = IV

switch flag
    case 1
        i=1;
        q0 = 0; % ug/hr
        V1x = 5000; % mL
        V2x = 5000; % mL
        CA0 = 13.7; % ug/ml (will be * ml = ug)
        kA =  .764; % hr-1

        OralOrIV = 1; % 1 = oral; 2 = IV
        [AUC1(i),AUC2(i),T1,Y1] = Acet_sim_t(q0,CA0,kA,V1x,V2x,OralOrIV,VisibleFigs);

        i=2;
        q0 = 0; % ug/hr
        V1x = 5000; % mL
        V2x = 5000; % mL
        CA0 = 12.9; % ug/ml * ml = ug
        kA =  1.27; % hr-1

        [AUC1(i),AUC2(i),T2,Y2] = Acet_sim_t(q0,CA0,kA,V1x,V2x,OralOrIV,VisibleFigs);
        
        AUC1
        AUC2
        
        figure;
        bar(AUC1);
        title(gca,'AUC for specific scenarios - central') 
        ylabel(gca,'AUC (ug*hr/ml)')
        xlabel(gca,'Scenario')
        set(gca,'XTickLabel',{'Softgel','Tablet'});
        
        figure;
        ax3=subplot(1,2,1);
        plot(ax3,T1,Y1(:,1),'k','linewidth',3);
        hold on;
        plot(ax3,T2,Y2(:,1),'r.','linewidth',3);
        title(ax3,'Concentration of Drug in Compartment 1') 
        ylabel(ax3,'[D] (ug/ml)')
        xlabel(ax3,'Time (hrs)')

        ax4=subplot(1,2,2);
        plot(ax4,T1,Y1(:,2),'k','linewidth',3);
        hold on;
        plot(ax4,T2,Y2(:,2),'r.','linewidth',3);
        title(ax4,'Concentration of Drug in Compartment 1') 
        ylabel(ax4,'[D] (ug/ml)')
        xlabel(ax4,'Time (hrs)')
        
    case 2

        i=1;
        q0 = 0; % ug/hr
        V1x = 5000; % mL
        V2x = 5000; % mL
        CA0 = 13.7; % ug/ml (will be * ml = ug)
        kA =  .764; % hr-1

        OralOrIV = 1; % 1 = oral; 2 = IV
        [AUC1(i),AUC2(i),T1,Y1] = Acet_sim_t(q0,CA0,kA,V1x,V2x,OralOrIV,VisibleFigs);

        i=2;
        OralOrIV = 2; % 1 = oral; 2 = IV
        [AUC1(i),AUC2(i),T2,Y2] = Acet_sim_t(q0,CA0,kA,V1x,V2x,OralOrIV,VisibleFigs);
        
        AUC1
        AUC2
        
        figure;
        bar(AUC1);
        title(gca,'AUC for specific scenarios - central') 
        ylabel(gca,'AUC (ug*hr/ml)')
        xlabel(gca,'Scenario')
        set(gca,'XTickLabel',{'Softgel','Softgel-equiv iv'});
        
        figure;
        ax3=subplot(1,2,1);
        plot(ax3,T1,Y1(:,1),'k','linewidth',3);
        hold on;
        plot(ax3,T2,Y2(:,1),'r.','linewidth',3);
        title(ax3,'Concentration of Drug in Compartment 1') 
        ylabel(ax3,'[D] (ug/ml)')
        xlabel(ax3,'Time (hrs)')

        ax4=subplot(1,2,2);
        plot(ax4,T1,Y1(:,2),'k','linewidth',3);
        hold on;
        plot(ax4,T2,Y2(:,2),'r.','linewidth',3);
        title(ax4,'Concentration of Drug in Compartment 1') 
        ylabel(ax4,'[D] (ug/ml)')
        xlabel(ax4,'Time (hrs)')
        
    case 3
        
        i=1;
        q0 = 0; % ug/hr
        V1x = 5000; % mL
        V2x = 5000; % mL
        CA0 = 12.9; % ug/ml * ml = ug
        kA =  1.27; % hr-1

        OralOrIV = 1; % 1 = oral; 2 = IV
        [AUC1(i),AUC2(i),T1,Y1] = Acet_sim_t(q0,CA0,kA,V1x,V2x,OralOrIV,VisibleFigs);

        i=2;
        
        OralOrIV = 2; % 1 = oral; 2 = IV
        [AUC1(i),AUC2(i),T2,Y2] = Acet_sim_t(q0,CA0,kA,V1x,V2x,OralOrIV,VisibleFigs);
        
        AUC1
        AUC2
        
        figure;
        bar(AUC1);
        title(gca,'AUC for specific scenarios - central') 
        ylabel(gca,'AUC (ug*hr/ml)')
        xlabel(gca,'Scenario')
        set(gca,'XTickLabel',{'Tablet','Tablet-equiv iv'});
        
        figure;
        ax3=subplot(1,2,1);
        plot(ax3,T1,Y1(:,1),'k','linewidth',3);
        hold on;
        plot(ax3,T2,Y2(:,1),'r.','linewidth',3);
        title(ax3,'Concentration of Drug in Compartment 1') 
        ylabel(ax3,'[D] (ug/ml)')
        xlabel(ax3,'Time (hrs)')

        ax4=subplot(1,2,2);
        plot(ax4,T1,Y1(:,2),'k','linewidth',3);
        hold on;
        plot(ax4,T2,Y2(:,2),'r.','linewidth',3);
        title(ax4,'Concentration of Drug in Compartment 1') 
        ylabel(ax4,'[D] (ug/ml)')
        xlabel(ax4,'Time (hrs)')
        

end
