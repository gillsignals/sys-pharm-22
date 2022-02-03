function dydt = m1C_eqns_Binding(t,y,p)

q=p(1);
V=p(2);
kc1=p(3);
kc2=p(4);
kc3=p(5);
kab=p(6);
kba=p(7);

dydt = zeros(3,1);    % make it a column vector (e.g. (3,1)

 dydt(1) = q/V - kc1*y(1) - kab*y(1)*y(2) + kba *y(3); %A
 dydt(2) =     - kc2*y(2) - kab*y(1)*y(2) + kba *y(3); %B
 dydt(3) =     - kc3*y(3) + kab*y(1)*y(2) - kba *y(3); %AB

