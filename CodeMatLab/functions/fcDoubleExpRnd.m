	function r = myDoubleExpRnd(mu,m,n); 

%	myDoubleExpRnd is Alan Barr's adaptation of the MatLab function exprnd 
%	to the case of a symmetric double-exponential distribution 
%	with mean of zero and mean Mu for positive values 
%	and -Mu for negative values. 
%
%	Syntax: r = myDoubleExpRnd(mu,m,n); 
%

%	Written 13 Feb 2009

	r=exprnd(mu,m,n); 
	iSign=rand(m,n); % uniformly distributed 0 to 1
	iNeg=find(iSign<0.5); % assumes no 0.5, which is VERY rare.
	r(iNeg)=-r(iNeg); 
	
