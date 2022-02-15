function dydt = CaffParax_eqns(t,y,p)
% One-compartment model of Caffeine, including absorption from gut and
% first-order clearance. Caffeine can also be metabolized to paraxanthine.

%% PARAMETERS 

q=p(1);  % infusion rate
V1=p(2); % central compartment volume
kc1=p(3);
kc2=p(4);
kcp=p(5);
ka=p(6);

%% EQUATIONS

dydt = zeros(5,1);    % make it a column vector (e.g. (3,1)
% 1 = caffeine in body (mg/L)
% 2 = paraxanthine in body (mg/L)
% 3 = cleared caffeine (mg) - note: amount not concentration
% 4 = cleared paraxanthine (mg) - note: amount not concentration
% 5 = caffeine in gut (mg) - note: amount not concentration

 dydt(1) = q/V1 + ka*y(5)/V1 - kc1*y(1) - kcp*y(1);
 dydt(2) =                   - kc2*y(2) + kcp*y(1);
 dydt(3) =                     kc1*y(1)*V1;
 dydt(4) =                     kc2*y(2)*V1;
 dydt(5) =      - ka*y(5);
