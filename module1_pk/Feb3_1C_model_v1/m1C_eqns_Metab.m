function dydt = m1C_eqns_Metab(t,y,p)
% This function defines equations to simulate a one-compartment model
% with a drug or prodrug administered by infusion that can undergo enzymatic 
% metabolism; both the drug and the metabolite can undergo clearance. In
% this code, we assume clearance and metabolism are first-order processes.
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
% solver to here; this will have three elements, one for each molecule)
% p = vector of parameters (we define this in our main code, and pass
% it to the ODE solver, which passes it to this function)

%% PARAMETERS
%
% The parameters needed for the equation are passed to this function in 
% the array p, so the next 5 lines are a way to assign those 5 values to
% named parameters and therefore make the equation(s) below easier to read
% and easier to check/debug. Take care with the order of the values in p.
q=p(1); % infusion rate
V=p(2); % volume
kc1=p(3); % clearance rate constant for molecule A
kc2=p(4); % clearance rate constant for molecule B
kab=p(5); % conversion rate constant for molecule A to molecule B

%% EQUATIONS
% initialize dydt to be the right size
% make it a column vector, i.e. (n,1) where n = number of equations
dydt = zeros(2,1);    % make it a column vector (e.g. (2,1))

% for metabolism, assume first order reaction 
%   (see low substrate concentration assumption of Michaelis-Menten)
 dydt(1) = q/V - kc1*y(1) - kab*y(1) ; % dA/dt
 dydt(2) =     - kc2*y(2) + kab*y(1) ; % dB/dt
