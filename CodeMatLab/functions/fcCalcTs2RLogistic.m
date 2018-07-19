function yhat = CalcTs2RLogistic(beta,x)

	b1=beta(1); b2=beta(2); b3=beta(3); yhat=b1./(1+exp(b2*(b3-x)));

