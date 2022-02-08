function dydt = m1C_eqns_Binding(t,y,p)
% This function defines equations to simulate a one-compartment model
% with a drug (A), a plasma protein to which the drug can bind (B), and the 
% drug-protein complex (AB). All three molecules can undergo clearance. 
% aside from drug infusion and the clearance processes, binding and unbinding
% are the other processes included in the equations.
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
%  solver to here; this will have three elements, one for each molecule)
% p = vector of parameters (we define this in our main code, and pass
%  it to the ODE solver, which passes it to this function)

%% PARAMETERS
%
% The parameters needed for the equation are passed to this function in 
% the array p, so the next 7 lines are a way to assign those 7 values to
% named parameters and therefore make the equation(s) below easier to read
% and easier to check/debug. Take care with the order of the values in p.
q=p(1); % infusion rate
V=p(2); % volume
kcA=p(3); % clearance rate constant for molecule A
kcB=p(4); % clearance rate constant for molecule B
kcAB=p(5); % clearance rate constant for molecule AB
kon=p(6);  % 'on rate' = forward rate constant for A binding B
koff=p(7); % 'off rate' = reverse rate constant for A unbinding B

%% EQUATIONS
% initialize dydt to be the right size
% make it a column vector, i.e. (n,1) where n = number of equations
dydt = zeros(3,1);  

% List of equations. Equation (1) is the ODE associated with the
% concentration of molecule (1); y(1) is the concentration of molecule (1).
% Note that the right-hand side of each equation depends on 
% the concentrations of multiple molecules. This is what makes these
% 'coupled ODEs'
%
% Also, note how appropriate spacing makes it easy 
% to see the symmetry across the various equations

dydt(1) = q/V - kcA*y(1)  - kon*y(1)*y(2) + koff*y(3); % dA/dt
dydt(2) =     - kcB*y(2)  - kon*y(1)*y(2) + koff*y(3); % dB/dt
dydt(3) =     - kcAB*y(3) + kon*y(1)*y(2) - koff*y(3); % dAB/dt

%% Alternate way of writing/organizing equations:
% In the equations above, we calculate the various rates in each equation, 
% even though some of the rates are the same. Instead, we can calculate 
% key repeated rates just once and then include them in the equations.
%
% There are some advantages to this, including the potential to make fewer
% mistakes through repeat typing, and making the equations slightly shorter
% and easier to read. However, they also may require careful reading
% depending on your choice of variable names or how many rate terms there are.
%
% You can comment out the three lines of equations above and 
%  uncomment the four equations below to try it out 
%  and check that you get the same answers:
%
%  r_binding = kon*y(1)*y(2) - koff*y(3);
%
%  dydt(1) = q/V - kcA*y(1)  - r_binding; % dA/dt
%  dydt(2) =     - kcB*y(2)  - r_binding; % dB/dt
%  dydt(3) =     - kcAB*y(3) + r_binding; % dAB/dt

% In some codes we may encounter, lower case 'v' is used instead of 'r' to
% denote rate; v stands for velocity, an alternate term for rate
