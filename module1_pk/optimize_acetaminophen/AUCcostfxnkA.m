function AUCcostout = AUCcostfxnkA(x1)

q0 = 0; % ug/hr
CA0 = 13.7; % ug/ml (will be * ml = ug)
kA =  x1; % hr-1
TimeLen = 8;

        [auc,t,y] = Acet_sim_opt(q0,CA0,kA,TimeLen);

        texp = [0.5,1,2,4,8];
        yexp = [5,6.5,5,3,1];

        for j=1:length(texp)
            teval = abs(t-texp(j));
            [tmin tindex] = min(teval);
            AUCcostout(j) = y(tindex) - yexp(j);
        end

end

