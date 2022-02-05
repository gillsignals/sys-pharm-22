function dydt = m1C_eqns(t,y,p)
% This function defines equations to simulate a 'toy model'
% (i.e. a model that doesn't represent a specific drug).
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
% the array p, so the next 3 lines are a way to assign those 3 values to
% named parameters and therefore make the equation(s) below easier to read
% and easier to check/debug. Take care with the order of the values in p.

q=p(1); % gather the parameters from p.
V=p(2); % there are easier ways to do this, we've written it this way 
k=p(3); % to emphasize the concept of passing parameters;
% later we'll demonstrate other ways to do this

%% EQUATIONS
% initialize dydt to be the right size
% make it a column vector, i.e. (n,1) where n = number of equations
dydt = zeros(1,1);    

% dydt(1) = q - k*y(1); % use this equation for q in units of nM/hr
dydt(1) = q/V - k*y(1) ; % use this equation for q in units of nmol/hr
 
