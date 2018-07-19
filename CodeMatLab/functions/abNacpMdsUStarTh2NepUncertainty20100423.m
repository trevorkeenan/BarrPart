function [Ensemble,xNEEgf] = ... 
	abNacpMdsUStarTh2NepUncertainty20100423 ... 
	(t,NEE,uStar,PPFD,Ta,Ts,Vpd,PPFDGF,TaGF,TsGF,VpdGF,fNight,uStarThs,iFig,cSiteYr); 

%abNacpMdsUStarTh2NepUncertainty20100423
%	estimates the uncertainty in gap-filled NEE 
%	that arises from uncertainty in the uStarTh. 
%
%Uncertainty in the uStarTh is estimated outside of this function 
%	by bootstrapping, and is input through the array: uStarThs. 
%	The number of uStarTh bootstrapping estimates (nBoot)
%	is typically set to 1,000. The uStarThs array may contain:
%-	single, annual uStarTh means (in a 1 x nBoot array), or
%-	annual sine curve fits, allowing for seasonal variation, 
%	(in a nt x nBoot array). 
%
%Syntax: 
%
%	[Ensemble,MatrixNEE] = abNacpMdsUStarTh2NepUncertainty20100423 ... 
%		(t,NEE,uStar,PPFD,Ta,Ts,PPFDGF,TaGF,TsGF,fNight,uStarThs);	
%
%	-	Ensemble is a structured record containing the ensemble output
%		(typically nBoot = 1,000 realizations) of NEE, RE and GPP, 
%		-	aggregated as annual, monthly, weekly and daily means, 
%		-	and averaged over the diel (24-h daily) cycle for 
%			annual, monthly and weekly periods. 
%	-	MatrixNEE contains all realizations (typically nBoot = 1,000) 
%		of the NEE time series. 
%
%	-	t is the time vector (with nt elements, usually one year)
%	-	NEE, uStar, PPFD, Ta, Ts and fNight are 
%		time series column vectors with nt elements:
%		-	NEE can be with or without u* filtering. 
%		-	PPFD is global incident radiation, either shortwave or PAR.
%		-	Ta and Ts are air and shallow soil temperature.
%		-	fNight is a time series vector with nt elements
%			that identifies daytime (0) or nighttime (1) periods. 
%	-	uStarThs contains nBoot realizations of the uStarTh filter, either: 
% 		-	single annual values in a 1 x nBoot array, or
% 		-	annual sine curve fits, allowing for seasonal variation, 
%			in a nt x nBoot array. 
%
%	-	iFig specifies the figure number for plotting (when iFig>0)
%	-	cSiteYr is a 10-character string containing the FLUXNET site code 
%		and year e.g. 'USHa1-2001'

%	========================================================================
%	========================================================================

%	Written 16 April 2010 by Alan Barr

%	=======================================================================
%	=======================================================================

	nt=length(t); [ntuStarTh,nBoot]=size(uStarThs);
	
	[y,m,d]=datevec(t); iYrs=fcx2rowvec(unique(y));
	nRecsPerDay=round(1/nanmedian(diff(t)));
	nMinsPerDay=24*60/nRecsPerDay; k=12*nMinsPerDay*60/1e6;
	
%	========================================================================
		
%	Process and fill gaps using bootstrapped ensemble of uStarThs. 
								
	FracSEBClosure=1; 
	
	MT=NaN*ones(nt,nBoot); xNEEgf=MT; 
	aMT=NaN*ones(nBoot,1); aNEP=aMT; aRE=aMT; aGPP=aMT; nNEP=aMT; 
		nMissNEP=aMT; nMissRE=aMT; nMissGPP=aMT;
		
%	Initialize averaging periods	
		
	nd=nt/nRecsPerDay; nw=52; nm=12;
	nMinsPerDay=24*60/nRecsPerDay;
	k=12*nMinsPerDay*60/1e6; % convert to g C m-2. 
	
	MatrixNEE=NaN*ones(nt,nBoot);
	
	aMT=NaN*ones(nBoot,1); aNEP=aMT; % a annual
	nMissNEP=aMT; % a annual
	
	dwMT=NaN*ones(nBoot,nRecsPerDay,nw); dwNEP=dwMT; % dw weekly mean diurnal cycle
	dmMT=NaN*ones(nBoot,nRecsPerDay,nm); dmNEP=dmMT; % dm monthly mean diurnal cycle
	daMT=NaN*ones(nBoot,nRecsPerDay); daNEP=daMT; % da annual mean diurnal cycle
	
	sdMT=NaN*ones(nBoot,nd); sdNEP=sdMT; % sd daily mean seasonal cycle
	swMT=NaN*ones(nBoot,nw); swNEP=swMT; % sw weekly mean seasonal cycle
	smMT=NaN*ones(nBoot,nm); smNEP=smMT; % sm monthly mean seasonal cycle
	
	clear *MT;
	
	itw=1:(nRecsPerDay*nw*7); % indices for weekly averaging, includes 52*7 days only. 
	
%	Estimate NEP at multiple uStarTh values.
	
	for iBoot=1:nBoot;
		
		uStarTh=uStarThs(:,iBoot); fPlot=0; 
		
		[NEEfilt,NEEgf,NEEhat] = ... 
			abNacpMdsCO2Flux2NEP20090205 ...
			(t,NEE,uStar,PPFDGF,TaGF,TsGF,VpdGF,PPFDGF,TaGF,TsGF,VpdGF,fNight,uStarTh,fPlot,cSiteYr);
		
		xNEEgf(:,iBoot)=NEEgf; NEPgf=-NEEgf; 
		
		%	Integrated annual seasonal cycles, g C m-2 t-1.
		
		aNEP(iBoot)=k*nansum(NEPgf); nNEP(iBoot)=sum(~isnan(NEPgf)); nMissNEP(iBoot)=sum(isnan(NEPgf));
		sdNEP(iBoot,:)=k*nansum(reshape(NEPgf,nRecsPerDay,nd));
		swNEP(iBoot,:)=k*nansum(reshape(NEPgf(itw),nRecsPerDay*7,nw));
		for im=1:12;
			it=(find(m==im)); nit=length(it);
			smNEP(iBoot,im)=k*nansum(NEPgf(it));
		end;
		
		%	Mean diurnal cycles, umol m-2 s-1.
		
		daNEP(iBoot,:,:)=nanmean(reshape(NEPgf,nRecsPerDay,nt/nRecsPerDay),2);
		dwNEP(iBoot,:,:)=nanmean(reshape(NEPgf(itw),nRecsPerDay,7,nw),2);
		for im=1:12;
			it=(find(m==im)); nit=length(it);
			dmNEP(iBoot,:,im)=nanmean(reshape(NEPgf(it),nRecsPerDay,nit/nRecsPerDay),2);
		end;

		disp(sprintf('Filling gaps in %s, %g/%g  %5.3f   %4.0f   n %5.0f nMiss %4.0f ', ...
			cSiteYr,iBoot,nBoot, mean(uStarTh), aNEP(iBoot), nNEP(iBoot), nMissNEP(iBoot)));
	
	end;
	
	Ensemble.NEE.nMissing = nMissNEP; 
	
	Ensemble.NEE.Annual = -aNEP; 
	Ensemble.NEE.Monthly = -smNEP; 
	Ensemble.NEE.Weekly = -swNEP; 
	Ensemble.NEE.Daily = -sdNEP; 
	
	Ensemble.NEE.DielByYr = -daNEP; 
	Ensemble.NEE.DielByMo = -dmNEP; 
	Ensemble.NEE.DielByWeek = -dwNEP; 
	
%	=======================================================================

	if iFig>0; 		
		
		fcFigLoc(iFig,0.3,0.45,'SW'); 
			
		a=-aNEP; cFlx=[cSiteYr ' NEE'];
		subplot('position',[0.15 0.76 0.82 0.18]); box on;
		hist(a,30); title(sprintf('%s %g (%g)',cFlx,round(mean(a)),round(std(a))));
		ylabel('Histogram (g C m^{-2})');
		
		sm=-smNEP; cFlx='NEE';
		sp=[0.15 0.41 0.82 0.28];
		subplot('position',sp); box on;
		boxplot(sm); set(gca,'position',sp);
		ylabel('Monthly (g C m^{-2})');
		
		da=-daNEP;
		sp=[0.15 0.06 0.82 0.28];
		subplot('position',sp); box on;
		boxplot(da); set(gca,'position',sp);
		ylabel('Diel (umol m^{-2} s^{-1})');
	
	end; 

%	=======================================================================
%	=======================================================================



