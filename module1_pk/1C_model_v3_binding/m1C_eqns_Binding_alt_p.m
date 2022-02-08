function dydt = m1C_eqns_Binding_alt_p(t,y,p)
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
% p = structured parameter set (we define this in our main code, and pass
%  it to the ODE solver, which passes it to this function)


%% Parameters
%
% This part is obsolete since we are importing p as a structure with named
% elements instead of as a vector of values
% q=p(1);
% V=p(2);
% kcA=p(3);
% kcB=p(4);
% kcAB=p(5);
% kon=p(6);
% koff=p(7);


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

 dydt(1) = p.q/p.V - p.kcA*y(1)  - p.kon*y(1)*y(2) + p.koff*y(3); %A
 dydt(2) =         - p.kcB*y(2)  - p.kon*y(1)*y(2) + p.koff*y(3); %B
 dydt(3) =         - p.kcAB*y(3) + p.kon*y(1)*y(2) - p.koff*y(3); %AB
 
 
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
%  r_binding = p.kon*y(1)*y(2) - p.koff *y(3);
%
%  dydt(1) = p.q/p.V - p.kcA*y(1)  - r_binding; %A
%  dydt(2) =         - p.kcB*y(2)  - r_binding; %B
%  dydt(3) =         - p.kcAB*y(3) + r_binding; %AB

% In some codes we may encounter, lower case 'v' is used instead of 'r' to
% denote rate; v stands for velocity, an alternate term for rate
 

