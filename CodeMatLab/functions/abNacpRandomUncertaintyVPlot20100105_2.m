	function [MuHat,b,cib,r2,s,mxNeg,myNeg,mxPos,myPos]= ... 
		abNacpRandomUncertaintyVPlot20100105_2 ...
		(NEE,NEEHat,iFig,cSiteYr,cFlux); 

    % this version outputs the parameters of the random uncertainty fit
    % lines displayed in the figure
    
%	function abNacpRandomUncertaintyVPlot20100105  
%
%	b=nacpRandomUncertaintyVPlot(NEE,NEEHat,iFig,cSiteYr); 
%
%	creates and fits the Andrew Richardson random uncertainty V plot 
%	based on measured and gap-filled modelled data. 
%
%	Version _20100105 uses regress to fit V plots 
%	and adds confidence intervals to the lin reg coefficients. 

%	Written 18 Nov 2009

%	========================================================================
%	========================================================================

	SqRt2=2^0.5; % convert double exponential Mu to Sigma equivalent
						
	Bias=NEEHat-NEE; iYaN=find(~isnan(Bias));
	
%	Compute bined means and fit Pos Neg separately
	
	nBins=30; ip=0:(100/nBins):100;
	
	MT=NaN*ones(nBins,1); mxNeg=MT; myNeg=MT; mxPos=MT; myPos=MT; 
	
	for iNP=1:2;
		
		switch iNP;
			case 1; it=find(NEEHat<0);
			case 2; it=find(NEEHat>=0);
		end;
		it=intersect(it,iYaN); 
		p=prctile(NEEHat(it),ip);
		for iBin=1:nBins;
			pL=p(iBin); pU=p(iBin+1);
			jt=find(NEEHat>=pL & NEEHat<=pU);
			jt=intersect(it,jt); njt=length(jt);
			Mu=expfit(abs(Bias(jt)));
			switch iNP;
				case 1; mxNeg(iBin)=mean(NEEHat(jt)); myNeg(iBin)=Mu;
				case 2; mxPos(iBin)=mean(NEEHat(jt)); myPos(iBin)=Mu;
			end;
		end; % for iBin
		
	end; % for iNP
	
%	Fit linear models. 
	
	[bNeg,bIntNeg,rNeg,rIntNeg,sNeg]=regress(myNeg,[ones(nBins,1) mxNeg]); 
		cibNeg=0.5*diff(bIntNeg,[],2); sNeg=std(rNeg); r2Neg=sNeg(1); 
	[bPos,bIntPos,rPos,rIntPos,sPos]=regress(myPos,[ones(nBins,1) mxPos]); 
		cibPos=0.5*diff(bIntPos,[],2); sPos=std(rPos); r2Pos=sPos(1); 
	
	b=[bNeg bPos]; cib=[cibNeg; cibPos]; r2=[r2Neg; r2Pos]; s=[sNeg; sPos]; 
	
% % % 	[mean(myNeg)/mean(mxNeg) mean(myPos)/mean(mxPos)]
		
	MuHat=NaN*NEEHat; 
	iNeg=find(NEEHat<0); MuHat(iNeg)=bNeg(1)+bNeg(2)*NEEHat(iNeg); 
	iPos=find(NEEHat>=0); MuHat(iPos)=bPos(1)+bPos(2)*NEEHat(iPos); 
	disp([fcx2rowvec(b) SqRt2*nanmean(MuHat)]);
	
	if iFig>0; 
	
		fcFigLoc(iFig,0.30,0.35,'NW'); 
		
		hold on;
		plot(mxNeg,SqRt2*myNeg,'ro');
		plot(mxPos,SqRt2*myPos,'bo');
		plot(NEEHat,SqRt2*MuHat,'k.','MarkerSize',3);
		hold off; grid on; box on;
		xlabel(cFlux); ylabel('2^{0.5}\mu');
		cTitle=sprintf('%s   %3.1f %4.2f   %3.1f %4.2f   %4.2f %4.2f ', ...
			cSiteYr,b,-mean(myNeg)/mean(mxNeg),mean(myPos)/mean(mxPos));
		title(cTitle);
		
	end; 
	
%	========================================================================
%	========================================================================


