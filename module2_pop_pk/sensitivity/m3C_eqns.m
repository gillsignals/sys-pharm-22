function dydt = m3C_eqns(t,y,p)

q=p(1);
V1=p(2);
V2=p(3);
V3=p(4);
kc1=p(5);
k12=p(6);
k21=p(7);
k13=p(8);
k31=p(9);
kDBon =p(10); % nM-1 hr-1
kDBoff=p(11);  %hr-1

dydt = zeros(6,1);    % make it a column vector (e.g. (3,1)

% 1 = drug in blood
% 2 = drug in tumor 1
% 3 = drug in body  
% 4 = drug in degradation compt (amount)
% 5 = binding protein in blood
% 6 = drug-binding protein; 

 dydt(1) = q/V1 - kc1*y(1)  - k12*y(1) + (V2/V1)*k21*y(2)  ...
                            - k13*y(1) + (V3/V1)*k31*y(3)  ...
                            - kDBon*y(1)*y(5) + kDBoff*y(6);
 dydt(2) =                  + (V1/V2)*k12*y(1) - k21*y(2);
 dydt(3) =                  + (V1/V3)*k13*y(1) - k31*y(3);
 dydt(4) =        kc1*y(1)*V1;
 dydt(5) =                  - kDBon*y(1)*y(5) + kDBoff*y(6) ;
 dydt(6) =                  + kDBon*y(1)*y(5) - kDBoff*y(6) ;
 
 
