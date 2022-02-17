close all;
clear all;

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

%% KEY SIMULATION PARAMETERS
flag=1;
% 1 = local univariate
% 2 = local bivariate
% 3 = local parameter vs global dose
% 4 = local parameter vs global plasma protein concentration

TimeLen=24;

OutputVar = 1;
    % 1 = drug in blood
    % 2 = drug in tumor
    % 3 = drug in body  
    % 4 = drug in degradation compt (amount)
    % 5 = binding protein in blood
    % 6 = drug-binding protein complex; 

ParamDelta = 0.05; % test sensitivity to a 5% change

%% PARAMETERS 

% Model parameters for array p
q = 0; % nmol/hr
V1ext = 0.5; % blood
V2ext = 0.2; % tumor
V3ext = 0.12; % rest
V1 = 5.0 * V1ext; % L
V2 = 0.5 * V2ext; % L
V3 = 40 * V3ext; % L
CL = .100; % L hr-1
kc1 = CL/V1; % hr-1

k12 = 0.1 * V2/V1;  % hr-1
k21 = 0.1; % hr-1
k13 = 0.5 * V3/V1;  % hr-1
k31 = 0.5 ; % hr-1

kDBoff = 0.1;  %hr-1
kDBon = kDBoff/6.5; % nM-1 hr-1   % Kd = 6.5 = koff/kon

% Initial concentrations
B0 = 50; %nM
D0 = 100; %nmol

%% DIVERGENT COLOR MAP
% Matlab doesn't have great divergent color maps (roughly speaking, spectra 
% with light colors in the middle and strong colors either side - so we
% add our own.
%
% This divergent colormap goes to red (up) and blue (down), with white in the middle
% May not be ideal or colorblind safe, hence the next version
% mapr = [linspace(0,1,101) ones(1,100)];
% mapg = [linspace(0,1,101) linspace(.99,0,100)];
% mapb = [ones(1,100) linspace(1,0,101)];
% map =[mapr' mapg' mapb'];
%
% This divergent color map goes to green (up) and purple (down), with white in the middle
% somewhat similar to examples on https://colorbrewer2.org/
mapr = [linspace(0.33,1,101) linspace(.99,0,100)]; 
mapg = [linspace(0,1,101) linspace(.99,.33,100)];
mapb = [linspace(0.33,1,101) linspace(0.99,0,100)];
map =[mapr' mapg' mapb'];


%% RUN BASE CASE

p0 = [q V1 V2 V3 kc1 k12 k21 k13 k31 kDBon kDBoff B0 D0]';
p0labels = {'q' 'V_1' 'V_2' 'V_3' 'k_c_1' 'k_1_2' 'k_2_1' 'k_1_3' 'k_3_1' 'k_D_B_o_n' 'k_D_B_o_f_f' 'B_0' 'D_0'}';

[auc0,t0,y0] = m3C_sim(p0,OutputVar,TimeLen,'on');
t = t0;
y = y0;

%% RUN SENSITIVITY SIMULATIONS

switch flag
    case 1
        
%% 1. SENSITIVITY - LOCAL UNIVARIATE
%========================
% OUTPUT: 1 = blood concentration of drug 
% OUTPUT: 6 = blood concentration of drug-protein complex
% INPUT: all parameters

for i=1:length(p0)
    p=p0;
    p(i)=p0(i)*(1.0+ParamDelta);
    [auc(i),t1,y1] = m3C_sim(p,OutputVar,TimeLen,'off');
    t = [t t1];
    y = [y y1];
    Sens(i) = ((auc(i)-auc0)/auc0)/((p(i)-p0(i))/p0(i)); % relative sensitivity vs parameter
    SensB(i) = ((auc(i)-auc0))/((p(i)-p0(i))); % absolute sensitivity vs parameter
    SensAbs(i) = ((auc(i)-auc0)/auc0); % relative change (not normalized to parameter)
    SensT = ((y1-y0)./y0)/((p(i)-p0(i))/p0(i)); % time-dependent relative sensitivity
    [maxS,indT] = max(abs(SensT)); 
    SensTmax(i) = t1(indT); % time to max relative sensitivity

    if i==1
        SensTarr = SensT;
    else
        SensTarr = [SensTarr SensT];
    end
end

figure;
hold on
for i=1:(length(p0)+1)
    plot(t(:,i),y(:,i));
end
ax = gca; % assign a handle (ax) for the current axes
ax.FontSize = 14;
ylabel('Concentration, Y')
xlabel('Time (hr)')
lgd = legend(p0labels);
lgd.Location = 'eastoutside';
lgd.Title.String = ['Parameter' newline 'being' newline 'changed'];


figure;
hold on
colors=lines(length(p0));
for i=1:length(p0)
    plot(t(:,i+1),SensTarr(:,i),'LineWidth',3);
    text(t(end,i+1)+.5,SensTarr(end,i),strcat(p0labels{i},{''}),'HorizontalAlignment','left','Color',colors(i,:));
end
ax = gca; % assign a handle (ax) for the current axes
ax.FontSize = 14;
ylabel({'Relative sensitivity of concentration','(dY/Y)/(dP/P)'})
xlabel('Time (hr)')
    
    figure;
    bar(Sens);  
	ax = gca; % assign a handle (ax) for the current axes
    ax.FontSize = 14;
    ylabel({'Relative sensitivity of AUC';'(dY/Y)/(dP/P)'})
    xlabel('Parameters')
    xticks([1:1:length(p0)])
    xticklabels(p0labels);
    xlim([.5 (length(p0)+.5)])

    figure;
    bar(SensB);
	ax = gca; % assign a handle (ax) for the current axes
    ax.FontSize = 14;
    ylabel({'Relative sensitivity of AUC';'(dY/dP)'})
    xlabel('Parameters')
    xticks([1:1:length(p0)])
    xticklabels(p0labels);
    xlim([.5 (length(p0)+.5)])
    
    figure;
    bar(SensTmax);
	ax = gca; % assign a handle (ax) for the current axes
    ax.FontSize = 14;
    ylabel({'Time of maximum sensitivity of concentration','(hr)'})
    xlabel('Parameters')
    xticks([1:1:length(p0)])
    xticklabels(p0labels);
    xlim([.5 (length(p0)+.5)])

    figure;
    subplot(2,1,1)
    bar(Sens);
	ax = gca; % assign a handle (ax) for the current axes
    ax.FontSize = 14;
    ylabel({'Relative sensitivity of AUC';'(dY/Y)/(dP/P)'})
    xlabel('Parameters')
    xticks([1:1:length(p0)])
    xticklabels(p0labels);
    xlim([.5 (length(p0)+.5)])
    subplot(2,1,2)
    bar(SensTmax);
	ax = gca; % assign a handle (ax) for the current axes
    ax.FontSize = 14;
    ylabel({'Time of maximum sensitivity','of concentration (hr)'})
    xlabel('Parameters')
    xticks([1:1:length(p0)])
    xticklabels(p0labels);
    xlim([.5 (length(p0)+.5)])
% figure;
% scatter([1:16],Sens,'k');


    case 2
%% 2. SENSITIVITY - LOCAL BIVARIATE
%========================
% OUTPUT: AUC of concentration
%         (if OutputVar = 1, then concentration in blood, etc)
% INPUT: all parameters x all parameters

for i=1:length(p0)
    p=p0;
    p(i)=p0(i)*(1.0+ParamDelta);
    [auc(i),t1,y1] = m3C_sim(p,OutputVar,TimeLen,'off');
    t = [t t1];
    y = [y y1];
    Sens(i) = ((auc(i)-auc0)/auc0)/((p(i)-p0(i))/p0(i));
    SensB(i) = ((auc(i)-auc0))/((p(i)-p0(i)));
    SensAbs(i) = ((auc(i)-auc0)/auc0);
    SensT = ((y1-y0)./y0)/((p(i)-p0(i))/p0(i));
end

for i=1:length(p0)
for j=1:length(p0)
    p=p0;
    p(i)=p0(i)*(1.0+ParamDelta);
    p(j)=p0(j)*(1.0+ParamDelta);
    [auc2(i,j),t1,y1] = m3C_sim(p,OutputVar,TimeLen,'off');
    Sens2(i,j) = ((auc2(i,j)-auc(j))/auc(j))/((p(i)-p0(i))/p0(i));
    Sens2Abs(i,j) = ((auc2(i,j)-auc0)/auc0);
    y0=y(:,j+1);
    Sens2T = ((y1-y0)./y0)/((p(i)-p0(i))/p0(i));
    [maxS,indT] = max(abs(Sens2T));
    Sens2Tmax(i,j) = t1(indT);
end
end

for i=1:length(p0)
for j=1:length(p0)
	Sens2Syn(i,j) = Sens2Abs(i,j)/(SensAbs(i)*SensAbs(j));
end
end

figure;
h1 = heatmap(p0labels,p0labels,Sens2);
h1.Colormap = map; % divergent colormap
h1.FontSize = 14;
h1.Title = 'Local sensitivity';
h1.XLabel = 'Parameter';
h1.YLabel = 'Parameter';

figure;
h2 = heatmap(p0labels,p0labels,Sens2Tmax,'CellLabelColor','none');
% h2.Colormap = map; % no need for divergent colormap
h2.FontSize = 14;
h2.Title = 'Time of peak sensitivity (hr)';
h2.XLabel = 'Parameter';
h2.YLabel = 'Parameter';

figure;
h3 = heatmap(p0labels,p0labels,Sens2Syn);
h3.Colormap = map; % divergent colormap
h3.FontSize = 14;
h3.Title = 'Synergy (Sij/(Si*Sj))';
h3.XLabel = 'Parameter';
h3.YLabel = 'Parameter';


    case 3
        
%% 3. SENSITIVITY - LOCAL BIVARIATE vs Dose
%========================
% OUTPUT: AUC of concentration
%         (if OutputVar = 1, then concentration in blood, etc)
% INPUT: all parameters and range of drug doses

tic 
p1 = p0;
Doses = round(logspace(1,4,10),0);
for i=1:length(Doses)
    % BASE CASE
p1(end) = Doses(i);
[auc0,t0,y0] = m3C_sim(p1,OutputVar,TimeLen,'off');

for j=1:length(p0)
    p=p1;
    p(j)=p1(j)*(1.0+ParamDelta);
    [auc3(i,j),t1,y1] = m3C_sim(p,OutputVar,TimeLen,'off');
    Sens3(i,j) = ((auc3(i,j)-auc0)/auc0)/((p(j)-p1(j))/p1(j));
    Sens3T = ((y1-y0)./y0)/((p(j)-p1(j))/p1(j));
    [maxS,indT] = max(abs(Sens3T));
    Sens3Tmax(i,j) = t1(indT);
end
end

toc
%Doselabels = num2str(Doses,'%6.0f');

tic
figure;
h4 = heatmap(Doses,p0labels,Sens3','CellLabelColor','none');
h4.Colormap = map;
h4.FontSize = 14;
h4.Title = 'Local (parameters) & Global (dose) sensitivity';
h4.XLabel = 'Dose (nmol)';
h4.YLabel = 'Parameter';

figure;
h5 = heatmap(Doses,p0labels,Sens3Tmax','CellLabelColor','none');
h5.FontSize = 14;
h5.Title = 'Time of peak sensitivity (hr)';
h5.XLabel = 'Dose (nmol)';
h5.YLabel = 'Parameter';
toc

    case 4

%% 4. SENSITIVITY - LOCAL BIVARIATE vs Plasma Protein level
%========================
% OUTPUT: AUC of concentration
%         (if OutputVar = 1, then concentration in blood, etc)
% INPUT: all parameters and range of drug doses

tic
p1 = p0;
Doses = round(logspace(1,4,10),0);
for i=1:length(Doses)
    % BASE CASE
p1(length(p1)-1) = Doses(i);
[auc0,t0,y0] = m3C_sim(p1,OutputVar,TimeLen,'off');

for j=1:length(p0)
    p=p1;
    p(j)=p1(j)*(1.0+ParamDelta);
    [auc3(i,j),t1,y1] = m3C_sim(p,OutputVar,TimeLen,'off');
    Sens3(i,j) = ((auc3(i,j)-auc0)/auc0)/((p(j)-p1(j))/p1(j));
    Sens3T = ((y1-y0)./y0)/((p(j)-p1(j))/p1(j));
    [maxS,indT] = max(abs(Sens3T));
    Sens3Tmax(i,j) = t1(indT);
end
end
toc

%Doselabels = num2str(Doses,'%6.0f');

tic
figure;
h6 = heatmap(Doses,p0labels,Sens3','CellLabelColor','none');
h6.Colormap = map;
h6.FontSize = 14;
h6.Title = 'Local (parameters) & Global (plasma protein) sensitivity';
h6.XLabel = 'Plasma Protein level (nM)';
h6.YLabel = 'Parameter';

figure;
h7 = heatmap(Doses,p0labels,Sens3Tmax','CellLabelColor','none');
h7.FontSize = 14;
h7.Title = 'Time of peak sensitivity (hr)';
h7.XLabel = 'Plasma Protein level (nM)';
h7.YLabel = 'Parameter';
toc 

end

% % Clustering - discuss later
% s3 = Sens3'
% s3(isnan(s3(:,1)),:)=[]
% x3 = pdist(s3);
% y3 = squareform(x3);
% z3 = linkage (y3);
% figure;
% dendrogram(z3);


