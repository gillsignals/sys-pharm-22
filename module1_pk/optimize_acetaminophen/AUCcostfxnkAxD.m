function AUCcostout = AUCcostfxnkAxD(x1)

q0 = 0; % ug/hr
CA0 = x1(2); % ug/ml (will be * ml = ug)
kA =  x1(1); % hr-1
TimeLen = 8;

        [auc,t,y] = Acet_sim_opt(q0,CA0,kA,TimeLen);

        texp = [0.5,1,2,4,8];
        yexp = [11,13,11,6,2];

        for j=1:length(texp)
            teval = abs(t-texp(j));
            [tmin tindex] = min(teval);
            AUCcostout(j) = y(tindex) - yexp(j);
        end

end

