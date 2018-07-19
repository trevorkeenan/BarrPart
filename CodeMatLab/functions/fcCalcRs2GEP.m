function yhat = CalcRs2GEP(beta,x)

GEPx=beta(1); bRs=beta(2);

Rs=x;


fRs=(Rs./(Rs+bRs));

yhat=GEPx*fRs;

if any(~isfinite(yhat))
    yhat=ones(length(x),1);
end




