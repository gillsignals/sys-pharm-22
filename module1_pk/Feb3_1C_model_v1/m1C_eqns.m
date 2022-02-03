function dydt = m1C_eqns(t,y,p)

%% Parameters 

q=p(1); % gather the parameters from p.
V=p(2); % there are easier ways to do this, we've written it this way 
k=p(3); % to emphasize the concept of passing parameters;
% later we'll demonstrate other ways to do this

dydt = zeros(1,1);    
% make it a column vector, i.e. (n,1) where n = number of equations

%% The Equations

% dydt(1) = q - k*y(1); % this equation for q in units of nM/hr
 dydt(1) = q/V - k*y(1) ; % this equation for q in units of nmol/hr
 
