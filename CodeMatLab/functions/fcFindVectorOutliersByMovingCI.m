	function [xClean,fOut]=FindVectorOutliersByMovingCI(t,x,nWindow,ns,fPlot); 
	
%	FindVectorOutliersByMovingCISameTimeOfDay
%
%	identifies outliers in a vector time series 

%	Updated 16 July 2009
	

	fCI='quartiles'; % options std iqr quartiles
	nHalfWindow=floor(nWindow/2); 
	fOut=zeros(size(x)); nOut=sum(fOut); nOut0=999; 
	
%	Allow for a number of iterations.

	xClean=x; nLoops=0; nLoopsX=1; 
	
	while nOut~=nOut0 && nLoops<nLoopsX; 
		
		nLoops=nLoops+1; nOut0=nOut; 
		iYaN=find(~isnan(xClean)); nYaN=length(iYaN); 
		
		for i=1:nYaN;
			i1=i-nHalfWindow; i2=i+nHalfWindow;
			if i1<1; i1=1; i2=nWindow; end;
			if i2>nYaN; i2=nYaN; i1=nYaN-nWindow; end;
               if i1<1; i1=1;end
			
				it=iYaN(i1:i2);
				switch fCI; % altered 6 July 2009.
					case 'std';
						m=mean(xClean(it)); s=std(xClean(it));
						if abs(xClean(iYaN(i))-m)>ns*s; fOut(iYaN(i))=1; end;
					case 'iqr';
						m=median(xClean(it)); s=iqr(xClean(it));
						if abs(xClean(iYaN(i))-m)>ns*s; fOut(iYaN(i))=1; end;
					case 'quartiles';
						q=prctile(xClean(it),[25 50 75]); 
						if (q(1)-x(iYaN(i)))>2*ns*(q(2)-q(1)); fOut(iYaN(i))=1; end;
						if (x(iYaN(i))-q(3))>2*ns*(q(3)-q(2)); fOut(iYaN(i))=1; end;
				end;
		end;
		iOut=find(fOut); nOut=length(iOut); xClean(iOut)=NaN;
		
	end;
	
	iOut=find(fOut);
	mx=nanmedian(x); iA=find(x(iOut)>mx); iB=find(x(iOut)<mx);
	nOut=length(iOut); nA=length(iA); nB=length(iB);
	disp(sprintf('Outliers:  %g %g %g ',[nOut nB nA])); 
	
	
	
	if fPlot>0;
		figure(fPlot); clf;
		subplot(2,1,1); plot(x,'.'); grid on; box on;
		f=fOut; f(find(f==0))=NaN; 
		hold on; plot(x.*f,'ko','MarkerFaceColor','c'); 
		subplot(2,1,2); plot(xClean,'.'); grid on; box on; 
		pause; 
	end; 
		
		
