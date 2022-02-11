function dydt = Acet_eqns(t,y,p)
% equations describing the pharmacokinetics of acetaminophen
% in a two-compartment model with absorption from the gut and
% clearance from the central compartment. IV infusion (q) is also included
% for comparison. A term for bioavailability is not included in these
% equations, but could be added.

dydt = zeros(4,1);    % use a column vector 
 
%% EQUATIONS
% 1 - concentration of drug in central compartment (infusion included here, 
%   but typically set to zero in simulations of oral delivery)
% 2 - concentration of drug in peripheral compartment (clearance included 
%   here, but typically set to zero in simulations)
% 3 - amount of drug in virtual clearance compartment. Note unit change
%   here, the terms in the equation have conc*vol form, which means amount 
%   not concentration. This eliminates the need to set a volume for this
%   compartment; just need to be careful in the main code to treat it like 
%   an amount, not a concentration.
% 4 - amount of drug in virtual gut compartment. Again, the units here are
%   in amount, not concentration. Note how y4 is divided by V in equation 1
%   to convert the amount back into a concentration in the central
%   compartment, using the central compartment volume as the correct
%   dilution

 dydt(1) = p.q/p.V1 + p.ka*y(4)/p.V1 - p.kc1*y(1) - p.k12*y(1) + (p.V2/p.V1)*p.k21*y(2);
 dydt(2) =                           - p.kc2*y(2) + (p.V1/p.V2)*p.k12*y(1) - p.k21*y(2);
 dydt(3) =                             p.kc1*y(1)*p.V1 + p.kc2*y(2)*p.V2;
 dydt(4) =          - p.ka*y(4);
 
 
