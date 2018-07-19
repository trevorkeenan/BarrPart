	function [eakPa,TdC]=myekPaTetenAboveWater(TaC,RHpc);

	eakPa=NaN*ones(size(TaC)); TdC=NaN*ones(size(TaC)); 
   
	a=1e-3*610.78; b=17.269; c=237.3; 
	eakPa=a*exp(b*TaC./(c+TaC)).*RHpc/100.; 
	d=log(eakPa/a)/b; TdC=d*c./(1-d); 

