function [out1,out2,out3] = Acet_sim_opt(q0,CA0,kA,TimeLen); 
% returns three outputs; receives four parameters. These can be modified as
% needed to return more/fewer outputs or receive more/fewer parameters. As 
% a general rule, pass the parameters/variables you think you will need,
% rather than passing everything.

% Acetaminophen - simulation code, for setting up the parameter values and
% calling the ode solver

%% PARAMETER VALUES
% note some of these values are passed in from the driver
% it's actually unnecessary to change the names (e.g. q = q0), you could
% just use the names above. The only reason I do it here is to emphasize
% where they come from, and to leave it flexible and easy to change which
% parameters get passed.

p.q = q0; % ug/hr
p.V1 = 5000; % mL
p.V2 = 5000; % mL
p.kc1 = .472; % hr-1
p.kc2 = 0; % hr-1
p.k12 = .319;  % hr-1
p.k21 = .499; % hr-1
p.ka =  kA; % hr-1

% Initial conditions
D0 = CA0*p.V1; % ug/ml * ml = ug
y0 = [0 0 0 D0]'; % y1,y2 in ug/ml; y3,y4 in ug

%% CALL SOLVER

options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
[T1,Y1] = ode45(@Acet_eqns,[0:TimeLen/1000:TimeLen],y0,options,p);

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
DrugIn = p.q*T1 + D0 ; % cumulative drug into system (continuous infusion + bolus)
DrugOut = CurrentDrug(:,3) ; % cumulative drug eliminated from system
BalanceD = DrugIn - DrugOut - CurrentDrug(:,1) - CurrentDrug(:,2) - CurrentDrug(:,4) + InitialDrug ; %(zero = balance)

% Add an automated check/report on mass balance, since we don't want to 
% look at hundreds of mass balance graphs
if mean(BalanceD)>1e-6
    disp('Mass imbalance possible: ');
    disp(BalanceD);
end

% calculate AUC by integrating the concentration curve (trapezoidal rule)
AUC = 0;
for i=1:(length(Y1)-1)
    AUC = AUC + 0.5*(Y1(i,1)+Y1(i+1,1))*(T1(i+1)-T1(i));
end

%% RETURN OUTPUTS

out1 = AUC;
out2 = T1;
out3 = Y1(:,1);

end
