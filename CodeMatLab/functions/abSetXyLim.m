	function fcSetXyLim(x,y,p,d); 
	
%	function fcSetXyLim(x,y,p,d)
%
%	sets x and y axis limits based on the p and 100-p percentiles, 
%	with d (fractional) spaces at the edge. 

	xN=prctile(x,p); xX=prctile(x,100-p); 
	yN=prctile(y,p); yX=prctile(y,100-p); 
	
	xx=x; iEx=find(xx<xN | xx>xX); xx(iEx)=NaN; 
	yy=y; iEx=find(yy<yN | yy>yX); yy(iEx)=NaN; 
	
	xN=min(xx); xX=max(xx); xR=xX-xN; 
	yN=min(yy); yX=max(yy); yR=yX-yN; 
	
	if xR>0; xlim([xN-d*xR xX+d*xR]); end; 
	if yR>0; ylim([yN-d*yR yX+d*yR]); end;
	


