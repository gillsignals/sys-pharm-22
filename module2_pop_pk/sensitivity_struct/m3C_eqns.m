function dydt = m3C_eqns(t,y,p)


dydt = zeros(6,1);    % make it a column vector (e.g. (3,1)

% 1 = drug in blood
% 2 = drug in tumor 1
% 3 = drug in body  
% 4 = drug in degradation compt (amount)
% 5 = binding protein in blood
% 6 = drug-binding protein; 

 dydt(1) = p.q/p.V1 - p.kc1*y(1)  - p.k12*y(1) + (p.V2/p.V1)*p.k21*y(2)  ...
                                  - p.k13*y(1) + (p.V3/p.V1)*p.k31*y(3)  ...
                            - p.kDBon*y(1)*y(5) + p.kDBoff*y(6);
 dydt(2) =                        + (p.V1/p.V2)*p.k12*y(1) - p.k21*y(2);
 dydt(3) =                        + (p.V1/p.V3)*p.k13*y(1) - p.k31*y(3);
 dydt(4) =            p.kc1*y(1)*p.V1;
 dydt(5) =                  - p.kDBon*y(1)*y(5) + p.kDBoff*y(6) ;
 dydt(6) =                  + p.kDBon*y(1)*y(5) - p.kDBoff*y(6) ;
 
 
