function yi=naninterp1(x,y,xi,Method);

	iYaN=find(~isnan(x+y)); nYaN=length(iYaN); 
	if nYaN<2; yi=NaN*ones(size(xi)); return; end; 
% % % 	[length(unique(x(iYaN))) length(unique(y(iYaN)))] 
% % % 	plot(x(iYaN),y(iYaN),'.'); 
   yi=interp1(x(iYaN),y(iYaN),xi,Method); nyi=length(yi);  
	
   iyiN=min(find(~isnan(yi))); if iyiN>1; yi(1:iyiN)=yi(iyiN); end; 
   iyiX=max(find(~isnan(yi))); if iyiX<nyi; yi(iyiX:nyi)=yi(iyiX); end; 
	
	iNaN=find(isnan(yi)); iLT=find(xi<x(1)); iLT=intersect(iLT,iNaN); yi(iLT)=y(1);  
	iNaN=find(isnan(yi)); iGT=find(xi>x(end)); iGT=intersect(iGT,iNaN); yi(iGT)=y(end);  
	

