function AUCcostout = AUCcostfxnD0(x1)

q0 = 0; % ug/hr
CA0 = x1; % ug/ml (will be * ml = ug)
kA =  .764; % hr-1
TimeLen = 8;

[auc,t,y] = Acet_sim_opt(q0,CA0,kA,TimeLen);

        texp = [0.5,1,2,4,8];
        yexp = [10,18,18,12,4];

        for j=1:length(texp)
            teval = abs(t-texp(j));
            [tmin tindex] = min(teval);
            AUCcostout(j) = y(tindex) - yexp(j);
        end

end

