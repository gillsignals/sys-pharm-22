function dydt = m1C_eqns_withDegr(t,y,p)
% This function defines equations to simulate a one-compartment model
% of one molecule - a drug delivered by infusion that can be cleared by a
% first-order process. In other words, the drug only undergoes two processes
% (infusion and clearance) and therefore only two rate terms appear in the 
% equation for the molecule. 
%
% Even though it's a one-compartment model, we add a second equation
% representing a second 'virtual' compartment, into which all the drug that
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
%  solver to here; this will have two elements, one for the drug in the 
%  blood (central compartment) and one for the cleared drug)
% p = structured parameter set (we define this in our main code, and pass
%  it to the ODE solver, which passes it to this function)

%% EQUATIONS
% initialize dydt to be the right size
% make it a column vector, i.e. (n,1) where n = number of equations
dydt = zeros(2,1);    

if p.cleared == 'conc'
  dydt(1) = p.q/p.V - p.k*y(1) ; % main compartment 
  dydt(2) =         + p.k*y(1) ; % degraded "compartment" 
 % Note that the drug in the virtual compartment (equation 2) is in units
 % of concentration. As a result, this compartment must have a volume
 % associated with it; we typically assume that it's the same volume as 
 % the central compartment
elseif p.cleared == 'amnt'
  dydt(1) = p.q/p.V - p.k*y(1) ; % main compartment 
  dydt(2) =         + p.k*y(1)*p.V ; % degraded "compartment" 
 % Note that the drug in the virtual compartment (equation 2) is in units
 % of amount - the rate of drug going into the compartment * volume. We
 % have to be careful with this, because the units are different in this
 % case - amount (e.g. mg, nmol) rather than concentration. We have to
 % account for that in the mass balance too.
else
    disp('error - no clearance compartment approach defined')
    return
end
