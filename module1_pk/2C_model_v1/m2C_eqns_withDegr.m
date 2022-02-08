function dydt = m2C_eqns_withDegr(t,y,p)
% This function defines equations to simulate a two-compartment model
% of one molecule - a drug delivered by infusion that can be cleared by a
% first-order process from the either compartment. The drug can also move 
% between the two compartments, again as a first-order process.
%
% Even though it's a two-compartment model, we add a third equation
% representing a third 'virtual' compartment, into which all the drug that
% is cleared/eliminated goes. By doing this, it makes it easy to keep a
% cumulative count of drug that has left the system (the compartment). 
%
% Our main code doesn't call this eqns function directly; 
% it calls an ODE solver, which then calls this function. 
%
%% THIS FUNCTION RETURNS:
% dydt = the rates of change of the concentration(s) of
%  the molecule(s) we are simulating, at time t
%
%% ARGUMENTS
% t = current time (this is passed from the ODE solver to here)
% y = current value of the concentrations (this is passed from the ode 
%  solver to here; this will have three elements, one for the drug in each
%  compartment and one for the cleared drug
% p = structured parameter set (we define this in our main code, and pass
%  it to the ODE solver, which passes it to this function)

%% EQUATIONS
% initialize dydt to be the right size
% make it a column vector, i.e. (n,1) where n = number of equations
dydt = zeros(3,1);    

%% EQUATIONS
%
% first attempt
%
dydt(1) = p.q/p.V1 - p.kc1*y(1) - p.k12*y(1) + p.k21*y(2);
dydt(2) =          - p.kc2*y(2) + p.k12*y(1) - p.k21*y(2);
dydt(3) =            p.kc1*y(1)*p.V1 + p.kc2*y(2)*p.V2; % NOTE: amount not concentration
 
%% EQUATIONS - updated with volume correction

% dydt(1) = p.q/p.V1 - p.kc1*y(1) - p.k12*y(1) + (p.V2/p.V1)*p.k21*y(2);
% dydt(2) =          - p.kc2*y(2) + (p.V1/p.V2)*p.k12*y(1) - p.k21*y(2);
% dydt(3) =            p.kc1*y(1)*p.V1 + p.kc2*y(2)*p.V2; % NOTE: amount not concentration
 
 
