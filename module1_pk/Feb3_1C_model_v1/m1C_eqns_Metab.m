function dydt = m1C_eqns_Metab(t,y,p)

q=p(1);
V=p(2);
kc1=p(3);
kc2=p(4);
kab=p(5);

dydt = zeros(2,1);    % make it a column vector (e.g. (3,1)
% assume first order reaction, low substrate
 dydt(1) = q/V - kc1*y(1) - kab*y(1) ; %A
 dydt(2) =     - kc2*y(2) + kab*y(1) ; %B
 
