function AUCcostout = AUCcostfxn1(x1)

q0 = 0; % ug/hr
CA0 = x1; % ug/ml (will be * ml = ug)
kA =  .764; % hr-1
TimeLen = 8;

[auc,t,y] = Acet_sim_opt(q0,CA0,kA,TimeLen);

AUCcostout = auc - 100;

end

