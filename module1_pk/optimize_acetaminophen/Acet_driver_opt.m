close all;
clear all;

RunCase = 1;

% 1. AUC depends on Dose
% 2. Optimization of dose to AUC
% 3. Optimization of dose to timepoint data
% 4. optimization with bounds (timepoint data)

% 5. 2d scan kcl x doses
% 6. Scan kcl vs one dose
% 7. optimization with bounds - absorption
% 8. optimization of two parameters - dose x absorption

q0 = 0; % ug/hr
CA0 = 13.7; % ug/ml (will be * ml = ug)
kA =  .764; % hr-1
TimeLen = 8;

% ParamDelta = 0.05; % test sensitivity to a 5% change

switch RunCase
    
    case 1
% 1. AUC depends on Dose
%========================

        Doses = CA0*logspace(-1,1,10);
        figure; 
        hold on;
        for i=1:length(Doses)
            [auc(i),t,y] = Acet_sim_opt(q0,Doses(i),kA,TimeLen);
            plot(t,y)
            title(gca,'Concentration of Drug in Blood')
            ylabel(gca,'[D] (ug/ml)')
            xlabel(gca,'time (hrs)')
        end
        lgd = legend(string(Doses));
        lgd.Location = 'best';
        lgd.Title.String = ['Doses ug/ml'];

        figure;
        hold on;
        scatter(Doses,auc)
        plot(Doses,auc)
        title(gca,'AUC vs Dose')
        xlabel(gca,'[D] (ug/ml)')
        ylabel(gca,'AUC (ug*hr/ml)')

    case 2
% 2. Optimization of dose to auc
%========================
     
        optimal = lsqnonlin(@AUCcostfxn1,5); % optimize
        % then run for these values!
        Dose = optimal;
        [auc,t,y] = Acet_sim_opt(q0,Dose,kA,TimeLen);
        auc
        Dose

    case 3
% 3. Optimization of dose to time points
%========================
     
        texp = [0.5,1,2,4,8];
        yexp = [10,18,18,12,4];
        optimal = lsqnonlin(@AUCcostfxnD0,10); % optimize
        % then run for these values!
        Dose = optimal;
        [auc,t,y] = Acet_sim_opt(q0,Dose,kA,TimeLen);
        auc
        Dose
        figure;
        hold on;
        plot(t,y)
        scatter(texp,yexp)
        title(gca,'Concentration of Drug in Blood')
        ylabel(gca,'[D] (ug/ml)')
        xlabel(gca,'time (hrs)')
    
    
	case 4
% 4. optimization with bounds 
%========================
     
        texp = [0.5,1,2,4,8];
        yexp = [10,18,18,12,4];
        x0=10;  % seed (initial guess)
        lb=1;   % lower bound
        ub=30;% upper bound
        options = optimoptions('lsqnonlin','Display','iter'); % display output
        optimal = lsqnonlin(@AUCcostfxnD0,x0,lb,ub,options)
        Dose = optimal;
        % check answer
        [auc,t,y] = Acet_sim_opt(q0,Dose,kA,TimeLen);
        
        figure;
        hold on;
        plot(t,y)
        scatter(texp,yexp)
        title(gca,'Concentration of Drug in Blood')
        ylabel(gca,'[D] (ug/ml)')
        xlabel(gca,'time (hrs)')
   

        
        
	case 5
% 5. 2d scan kA vs doses
%========================
     
        kAs = kA*logspace(-2,2,10); % hr-1

        for j=1:length(kAs)
            Doses = CA0*logspace(-1,1,10);
            for i=1:length(Doses)
                [auc(i,j),t,y] = Acet_sim_opt(q0,Doses(i),kAs(j),TimeLen);
            end
        end
        
        figure;
        hold on
        for j=1:length(kAs)
            plot(Doses,auc(:,j))
        end
        title(gca,'AUC vs Dose')
        xlabel(gca,'[D] (ug/ml)')
        ylabel(gca,'AUC (ug*hr/ml)')
        lgd = legend(string(kAs));
        lgd.Location = 'best';
        lgd.Title.String = ['kA values'];

        
        
	case 6
% 6. Scan kA vs one dose
%========================
        
        kAs = kA*logspace(-2,2,10); % hr-1

        figure;
        hold on;

        for j=1:length(kAs)
            [auc(j),t,y] = Acet_sim_opt(q0,CA0,kAs(j),TimeLen);
            plot(t,y)
        end
        title(gca,'Concentration of Drug in Blood')
        ylabel(gca,'[D] (ug/ml)')
        xlabel(gca,'time (hrs)')
        lgd = legend(string(kAs));
        lgd.Location = 'best';
        lgd.Title.String = ['kA values'];

        figure;
        subplot(1,2,1)
        hold on;
        scatter(kAs,auc)
        plot(kAs,auc)
        set(gca,'xscale','log')
        title(gca,'AUC vs absorption (log)')
        xlabel(gca,'k_A (hr^-^1)')
        ylabel(gca,'AUC (ug*hr/ml)')

        subplot(1,2,2)
        hold on;
        scatter(kAs,auc)
        plot(kAs,auc)
        set(gca,'xscale','linear')
        title(gca,'AUC vs absorption (linear)')
        xlabel(gca,'k_A (hr^-^1)')
        ylabel(gca,'AUC (ug*hr/ml)')


	case 7
% 7. optimization with bounds - kA
%========================

        texp = [0.5,1,2,4,8];
        yexp = [5,6.5,5,3,1];
        x0=kA;
        lb=kA/10;
        ub=kA*10;
        options = optimoptions('lsqnonlin','Display','iter');
        optimal = lsqnonlin(@AUCcostfxnkA,x0,lb,ub,options)
        % check answer
        [auc,t,y] = Acet_sim_opt(q0,CA0,optimal,TimeLen);
        
        figure;
        hold on;
        plot(t,y)
        scatter(texp,yexp)
        title(gca,'Concentration of Drug in Blood')
        ylabel(gca,'[D] (ug/ml)')
        xlabel(gca,'time (hrs)')

	case 8
% 8. optimization of two parameters - kA, Dose
%========================

        texp = [0.5,1,2,4,8];
        yexp = [11,13,11,6,2];
        x0=[kA CA0];
        lb=x0/10;
        ub=x0*10;
        options = optimoptions('lsqnonlin','Display','iter');
        optimal = lsqnonlin(@AUCcostfxnkAxD,x0,lb,ub,options)
        % check answer
        [auc,t,y] = Acet_sim_opt(q0,optimal(2),optimal(1),TimeLen);
        
        figure;
        hold on;
        plot(t,y)
        scatter(texp,yexp)
        title(gca,'Concentration of Drug in Blood')
        ylabel(gca,'[D] (ug/ml)')
        xlabel(gca,'time (hrs)')   
           
end
