function [Ensemble,xNEEgf] = ... 
	abNacpFcrnUStarTh2NepUncertainty20100423 ... 
	(t,NEE,uStar,PPFD,Ta,Ts,PPFDGF,TaGF,TsGF,fNight,uStarThs,cSiteYr,iFig,k); 

%abNacpFcrnUStarTh2NepUncertainty20100423
%	estimates the uncertainty in gap-filled NEE, RE and GPP 
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
%	[Ensemble,MatrixNEE] = abNacpFcrnUStarTh2NepUncertainty20100423 ... 
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

% 		try;
			
	nt=length(t); [ntuStarTh,nBoot]=size(uStarThs);
	
	if ~ismember(ntuStarTh,[nt 1]);
% % % 		assign outputs
% % % 		err msg
% % % 		return;
	end;
	
	
	[y,m,d]=datevec(t); iYrs=fcx2rowvec(unique(y));
	nRecsPerDay=round(1/nanmedian(diff(t)));
    
    % this section was giving errors from using round, tjk August 2012
    if nt<=8785    % to account for a leap year
        nRecsPerDay=24;
    else
        nRecsPerDay=48;
    end
    
	nMinsPerDay=24*60/nRecsPerDay; %k=12*nMinsPerDay*60/1e6;
	
%	========================================================================
		
%	Process and fill gaps using bootstrapped ensemble of uStarThs. 
								
	FracSEBClosure=1; 
	
	MT=NaN*ones(nt,nBoot); xNEEgf=MT; 
	aMT=NaN*ones(nBoot,1); aNEP=aMT; aRE=aMT; aGPP=aMT; nNEP=aMT; 
		nMissNEP=aMT; nMissRE=aMT; nMissGPP=aMT;
		
% % % 	add other averaging periods	
		
	nd=nt/nRecsPerDay; nw=52; nm=12;
	nMinsPerDay=24*60/nRecsPerDay;
	k=12*nMinsPerDay*60/1e6; % convert to g C m-2. 
	
	MatrixNEE=NaN*ones(nt,nBoot);
	
	aMT=NaN*ones(nBoot,1); aNEP=aMT; aRE=aMT; aGPP=aMT; % a annual
	nMissNEP=aMT; nMissRE=aMT; nMissGPP=aMT; % a annual
	
	dwMT=NaN*ones(nBoot,nRecsPerDay,nw); dwNEP=dwMT; dwRE=dwMT; dwGPP=dwMT; % dw weekly mean diurnal cycle
	dmMT=NaN*ones(nBoot,nRecsPerDay,nm); dmNEP=dmMT; dmRE=dmMT; dmGPP=dmMT; % dm monthly mean diurnal cycle
	daMT=NaN*ones(nBoot,nRecsPerDay); daNEP=daMT; daRE=daMT; daGPP=daMT; % da annual mean diurnal cycle
	
	sdMT=NaN*ones(nBoot,nd); sdNEP=sdMT; sdRE=sdMT; sdGPP=sdMT; % sd daily mean seasonal cycle
	swMT=NaN*ones(nBoot,nw); swNEP=swMT; swRE=swMT; swGPP=swMT; % sw weekly mean seasonal cycle
	smMT=NaN*ones(nBoot,nm); smNEP=smMT; smRE=smMT; smGPP=smMT; % sm monthly mean seasonal cycle
	
	clear *MT;
	
	itw=1:(nRecsPerDay*nw*7); % indices for weekly averaging, includes 52*7 days only. 
	
%	Estimate NEP at multiple uStarTh values.
	
	for iBoot=1:nBoot;
		
% % % 		sum(isnan(NEE))
		
		uStarTh=uStarThs(:,iBoot); fPlot=0; 
		
		[NEP,RE,GPP,NEPgf,REgf,GPPgf,RHat,GPPHat,RHat0,GPPHat0,mtR,mcR,mtGPP,mcGPP,bR,bGPP] = ... 
			abNacpFcrnCO2Flux2NEP20090205 ...
			(t,NEE,uStar,PPFDGF,TaGF,TsGF,PPFDGF,TaGF,TsGF,fNight,uStarTh,FracSEBClosure,cSiteYr,fPlot);
		
		xNEEgf(:,iBoot)=-NEPgf; 
		
		%	Integrated annual seasonal cycles, g C m-2 t-1.
		
		aNEP(iBoot)=k*nansum(NEPgf); nNEP(iBoot)=sum(~isnan(NEPgf)); nMissNEP(iBoot)=sum(isnan(NEPgf));
		aRE(iBoot)=k*nansum(REgf); nMissRE(iBoot)=sum(isnan(REgf));
		aGPP(iBoot)=k*nansum(GPPgf); nMissGPP(iBoot)=sum(isnan(GPPgf));
		
		sdNEP(iBoot,:)=k*nansum(reshape(NEPgf,nRecsPerDay,nd));
		sdRE(iBoot,:)=k*nansum(reshape(REgf,nRecsPerDay,nd));
		sdGPP(iBoot,:)=k*nansum(reshape(GPPgf,nRecsPerDay,nd));
		
		swNEP(iBoot,:)=k*nansum(reshape(NEPgf(itw),nRecsPerDay*7,nw));
		swRE(iBoot,:)=k*nansum(reshape(REgf(itw),nRecsPerDay*7,nw));
		swGPP(iBoot,:)=k*nansum(reshape(GPPgf(itw),nRecsPerDay*7,nw));
		
		for im=1:12;
			it=(find(m==im)); nit=length(it);
			smNEP(iBoot,im)=k*nansum(NEPgf(it));
			smRE(iBoot,im)=k*nansum(REgf(it));
			smGPP(iBoot,im)=k*nansum(GPPgf(it));
		end;
		
		%	Mean diurnal cycles, umol m-2 s-1.
		
		daNEP(iBoot,:,:)=nanmean(reshape(NEPgf,nRecsPerDay,nt/nRecsPerDay),2);
		daRE(iBoot,:,:)=nanmean(reshape(REgf,nRecsPerDay,nt/nRecsPerDay),2);
		daGPP(iBoot,:,:)=nanmean(reshape(GPPgf,nRecsPerDay,nt/nRecsPerDay),2);
		
		dwNEP(iBoot,:,:)=nanmean(reshape(NEPgf(itw),nRecsPerDay,7,nw),2);
		dwRE(iBoot,:,:)=nanmean(reshape(REgf(itw),nRecsPerDay,7,nw),2);
		dwGPP(iBoot,:,:)=nanmean(reshape(GPPgf(itw),nRecsPerDay,7,nw),2);
		
		for im=1:12;
			it=(find(m==im)); nit=length(it);
			dmNEP(iBoot,:,im)=nanmean(reshape(NEPgf(it),nRecsPerDay,nit/nRecsPerDay),2);
			dmRE(iBoot,:,im)=nanmean(reshape(REgf(it),nRecsPerDay,nit/nRecsPerDay),2);
			dmGPP(iBoot,:,im)=nanmean(reshape(GPPgf(it),nRecsPerDay,nit/nRecsPerDay),2);
		end;

		disp(sprintf('Filling gaps in %s, %g/%g  %5.3f   %4.0f %4.0f %4.0f   n %5.0f nMiss %4.0f %4.0f %4.0f', ...
			cSiteYr,iBoot,nBoot, mean(uStarTh), aNEP(iBoot),aRE(iBoot),aGPP(iBoot), nNEP(iBoot), nMissNEP(iBoot),nMissRE(iBoot),nMissGPP(iBoot)));
	
	end;
	
	Ensemble.NEE.nMissing = nMissNEP; Ensemble.RE.nMissing = nMissRE; Ensemble.GPP.nMissing = nMissGPP; 
	
	Ensemble.NEE.Annual = -aNEP; Ensemble.RE.Annual = aRE; Ensemble.GPP.Annual = aGPP; 
	Ensemble.NEE.Monthly = -smNEP; Ensemble.RE.Monthly = smRE; Ensemble.GPP.Monthly = smGPP; 
	Ensemble.NEE.Weekly = -swNEP; Ensemble.RE.Weekly = swRE; Ensemble.GPP.Weekly = swGPP; 
	Ensemble.NEE.Daily = -sdNEP; Ensemble.RE.Daily = sdRE; Ensemble.GPP.Daily = sdGPP; 
	
	Ensemble.NEE.DielByYr = -daNEP; Ensemble.RE.DielByYr = daRE; Ensemble.GPP.DielByYr =  daGPP; 
	Ensemble.NEE.DielByMo = -dmNEP; Ensemble.RE.DielByMo = dmRE; Ensemble.GPP.DielByMo = dmGPP; 
	Ensemble.NEE.DielByWeek = -dwNEP; Ensemble.RE.DielByWeek = dwRE; Ensemble.GPP.DielByWeek = dwGPP; 
	
%	=======================================================================

	if iFig>0; 		
		
		fcFigLoc(iFig,0.5,0.45,'SW'); 
		
		for iFlx=1:3; 
			
			switch iFlx; case 1; a=-aNEP; cFlx=[cSiteYr ' NEE']; case 2; a=aRE; cFlx='RE'; case 3; a=aGPP; cFlx='GPP'; end; 
			subplot('position',[0.08+(iFlx-1)*0.32 0.76 0.26 0.18]); box on; 
			hist(a,30); title(sprintf('%s %g (%g)',cFlx,round(mean(a)),round(std(a)))); 
			if iFlx==1; ylabel('Histogram (g C m^{-2})'); end; 
		
			switch iFlx; case 1; sm=-smNEP; cFlx='NEE'; case 2; sm=smRE; cFlx='RE'; case 3; sm=smGPP; cFlx='GPP'; end; 
			sp=[0.08+(iFlx-1)*0.32 0.41 0.26 0.28]; 
			subplot('position',sp); box on; 
			boxplot(sm); set(gca,'position',sp); 
			if iFlx==1; ylabel('Monthly (g C m^{-2})'); end; 
			
			switch iFlx; case 1; da=-daNEP; case 2; da=daRE; case 3; da=daGPP; end; 
			sp=[0.08+(iFlx-1)*0.32 0.06 0.26 0.28]; 
			subplot('position',sp); box on; 
			boxplot(da); set(gca,'position',sp); 
			if iFlx==1; ylabel('Diel (umol m^{-2} s^{-1})'); end; 
	
		end; 	

	end; 

%	=======================================================================
%	=======================================================================



