	function [xClean,fOut]=FindVectorOutliersByMovingPapaleSpikes(t,x,nWindow,ns,fPlot); 
	
%	FindVectorOutliersByMovingPapaleSpikes
%
%	identifies outliers 

%	Written 16 July 2009

	fCI='std'; % options std iqr quartiles

	x=fcx2colvec(x); xClean=x; xDirty=x; nt=length(t); fOut=zeros(size(xClean)); 
	
	nHalfWindow=floor(nWindow/2); nOut=sum(fOut); nOut0=999; 
	nLoops=0; nLoopsX=1; 
	
	while nOut~=nOut0 && nLoops<nLoopsX; 
		
		nLoops=nLoops+1; nOut0=nOut; 

		x0=[xClean(2:end); NaN];
		x2=[NaN; xClean(1:(end-1))];
		d=xClean-(x0+x2)/2;
		
		iYaN=find(~isnan(d)); nYaN=length(iYaN); 

		for i=1:nYaN;
			i1=i-nHalfWindow; i2=i+nHalfWindow;
			if i1<1; i1=1; i2=nWindow; end;
			if i2>nYaN; i2=nYaN; i1=nYaN-nWindow; end;
				it=iYaN(i); itWindow=iYaN(i1:i2);
				switch fCI; % altered 6 July 2009.
					case 'std';
						md=mean(d(itWindow)); sd=std(d(itWindow)); % disp([md sd]); 
						if abs(d(it)-md)>ns*sd; fOut(it)=1; end;
					case 'iqr';
						md=median(d(itWindow)); sd=iqr(d(itWindow));
						if abs(d(it)-md)>ns*sd; fOut(it)=1; end;
					case 'quartiles';
						q=prctile(d(itWindow),[25 50 75]); 
						if (q(1)-d(it))>2*ns*(q(2)-q(1)); fOut(it)=1; end;
						if (d(it)-q(3))>2*ns*(q(3)-q(2)); fOut(it)=1; end;
				end;
		end;
		iOut=find(fOut); nOut=length(iOut); xClean(iOut)=NaN;
		
	end;
		
	iOut=find(fOut);
	mx=nanmedian(x); iA=find(x(iOut)>mx); iB=find(x(iOut)<mx);
	nOut=length(iOut); nA=length(iA); nB=length(iB);
	disp(sprintf('Papale spikes:  %g %g %g ',[nOut nB nA])); 

	if fPlot; 
		figure(1); clf; 
		subplot(2,1,1); plot(xDirty,'.'); grid on; box on; 
		f=fOut; f(find(f==0))=NaN; 
		hold on; plot(xDirty.*f,'ko','MarkerFaceColor','y'); 
		subplot(2,1,2); plot(xClean,'.'); grid on; box on; 
		pause; 
	end; 
		
	