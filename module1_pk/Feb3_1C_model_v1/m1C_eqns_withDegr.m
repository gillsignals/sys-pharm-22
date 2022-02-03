function dydt = m1C_eqns_withDegr(t,y,p)

q=p(1);
V=p(2);
k=p(3);

dydt = zeros(2,1);    % make it a column vector (e.g. (3,1)

 dydt(1) = q/V - k*y(1) ;
 dydt(2) =     + k*y(1) ;
 
