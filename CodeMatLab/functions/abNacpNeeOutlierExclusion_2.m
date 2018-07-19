	function [NeeQC,iOut] = abNacpNeeOutlierExclusion_2 ... 
		(t,NEE,PAR,Ta,Ts,fNight,nsP,nsCi,ndWindow,cSiteYr,iFig,cFlux); 

%abNacpNeeOutlierExclusion
%	finds and excludes outliers in annual NEE time series 
%	using several approaches:
%	-	Papale (2006) spike detection
%	-	Papale (2006) spike detection modified to handle cases 
%		where adjacent points are missing (using modelled NEE). 
%	-	confidence intervals applied to the NEE gap-filling model 
%		residuals, implemented using moving windows to deal with 
%		seasonal changes in variance. 
%
%Syntax:
%
%	[NeeQC,iOut] = barrNacpNeeOutlierExclusion ... 
% 		(t,NEE,uStar,PAR,Ta,T,fNight,ns,ndWindow); 
%
%	-	NeeQC sets outliers to NaN
%	-	iOut is an index of excluded values
%
%	-	t is the time vector (with nt elements)
%	-	NEE, PAR, Ta, Ts and fNight are 
%		time series column vectors with nt elements:
%		-	NEE can be with or without u* filtering. 
%		-	PAR is global incident radiation, either shortwave or PAR.
%		-	Ta and Ts are air and shallow soil temperature.
%		-	fNight is a time series vector with nt elements
%			that identifies daytime (0) or nighttime (1) periods. 
%	-	uStarTh is the uStarTh filter, either: 
% 		-	a scalar (single annual value) 
% 		-	a vector with nt elements (allowing seasonal variation)
%		- typically set to zero for outlier selection. 
%
%	-	ns controls the strength of the outlier exclusion 
%		e.g., ns=7 (recommended) is mild, 5 is moderate, and 3 is strong. 
%	-	ndWindow controls the width (in days) of the moving windows 
%		in the analysis of confidence intervals.

%	=======================================================================
%	=======================================================================
	
	nRecsPerDay=round(1/nanmedian(diff(t))); 
	ntWindow=nRecsPerDay*ndWindow; 
	
	uStar=zeros(size(t)); uStarTh=0; FracSEBClosure=1; fPlot=0; 
	
	NeeQC=NEE; 
	
	% pass 1 on NEE, spikes only
	
	[junk,fOut]=fcFindVectorOutliersByMovingPapaleSpikes(t,NEE,ntWindow,nsP,fPlot);
	iOutNeeP1=find(fOut); NeeQC(iOutNeeP1)=NaN;
	
	% pass 2 on ResidNEP, spikes and confidence interval outliers
	
	[NEP,R,GPP,NEPgf,Rgf,GPPgf,RHat,GPPHat] = abNacpFcrnCO2Flux2NEP20090205 ...
		(t,NeeQC,uStar,PAR,Ta,Ts,PAR,Ta,Ts,fNight,uStarTh,FracSEBClosure,cSiteYr,fPlot);
	NEPHat=GPPHat-RHat; ResidNEP=(-NeeQC)-NEPHat;
	
	[junk,fOut]=fcFindVectorOutliersByMovingPapaleSpikesWrtxHat(t,NeeQC,-NEPHat,ntWindow,nsP,fPlot);
		iOutNeeP2=find(fOut); NeeQC(iOutNeeP2)=NaN; ResidNEP=(-NeeQC)-NEPHat;
	
	gfResidNEP=ResidNEP; gfResidNEP(isnan(gfResidNEP))=0; 	
	[junk,fOutResidP]=fcFindVectorOutliersByMovingPapaleSpikes(t,gfResidNEP,ntWindow,nsP,fPlot);
		iOutResidP=find(fOutResidP); NeeQC(iOutResidP)=NaN; ResidNEP=(-NeeQC)-NEPHat;
	
	[junk,fOutResidCi]=fcFindVectorOutliersByMovingCI(t,ResidNEP,ntWindow,nsCi,fPlot);
		iOutResidCi=find(fOutResidCi); NeeQC(iOutResidCi)=NaN;
	
	iOutP=union(iOutNeeP1,iOutNeeP2); 	
	iOut=union(iOutP,iOutResidCi);
	iOut=union(iOut,iOutResidP); 
	
	if iFig>0; 
		fcFigLoc(iFig,0.3,0.35,'NE'); MS='MarkerSize'; MFC='MarkerFaceColor'; 
		plot(t,NEE,'b.', t(iOut),NEE(iOut),'r.'); 
		fcDateTick(t,'Mo',4,1); title([cFlux ' Outliers ' cSiteYr]); 
	end; 
	
%	=======================================================================
%	=======================================================================
	
	