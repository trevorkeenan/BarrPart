	function [Ensemble,MatrixNEE,bNeg,bPos]=abNacpFcrnRandomUncertainty20100423_FLUX2 ...
		(t,NEE,NEEhat,uStar,Rsd,Ta,Ts,fNight,uStarTh,DirOut,cSiteYr,nMC,iFig,k,cFlux) 

    % this version outputs the parameters of the random uncertainty for
    % positive and negative bins

%abNacpFcrnRandomUncertainty20100423
%	random uncertainty analysis function, 
%	based on Richardson & Hollinger (2007), 
%	using a modified version of the Fluxnet-Canada 
%	gap-filling procedure (Barr et al. 2004), 
%	typically run one year at a time. 
%See also adrRandomUncertaintyMDS20100423.	
%
%Syntax: 
%
%	[Ensemble,MatrixNEE] = abNacpFcrnRandomUncertainty20100423 ...
% 		(t,NEE,uStar,Rsd,Ta,Ts,fNight,uStarTh,DirOut,cSiteYr,nMC); 
%
%	-	Ensemble is a structured record containing the ensemble output
%		(typically nMC = 1,000 realizations) of NEE, RE and GPP, 
%		-	aggregated as annual, monthly, weekly and daily means, 
%		-	and averaged over the diel (24-h daily) cycle for 
%			annual, monthly and weekly periods. 
%	-	MatrixNEE contains all realizations (typically nMC = 1,000) 
%		of the NEE time series. 
%
%	-	t is the time vector (with nt elements)
%	-	NEE, uStar, Rsd, Ta, Ts and fNight are 
%		time series column vectors with nt elements:
%		-	NEE can be with or without u* filtering. 
%		-	Rsd is global incident radiation, either shortwave or PAR.
%		-	Ta and Ts are air and shallow soil temperature.
%		-	fNight is a time series vector with nt elements
%			that identifies daytime (0) or nighttime (1) periods. 
%	-	uStarTh is the uStarTh filter, either: 
% 		-	a scalar (single annual value) 
% 		-	a vector with nt elements (allowing seasonal variation)
%
%	-	DirOut is the root directory for data input and output.
%	-	cSiteYr is a 10-character string containing the FLUXNET site code 
%		and year, separated by a dash e.g. 'USHa1-2000'
%	-	nMC is the number of Monte-Carlo realizations, typically 1,000. 
%

%	========================================================================
%	========================================================================

%	Functions called: 
%
%		fcBin, fcDateTick, fcEqnAnnualSine, fcNaniqr, fcReadFields 
%		fcr2Calc, fcx2colvec, fcx2rowvec
%		stats toolbox:  nanmedian, nanmean, nlinfit, prctile

%	fcDateVec

%	Written 16 April 2010 by Alan Barr

%	=======================================================================
%	=======================================================================

%	Estimate NEEhat (modelled NEE from Fluxnet-Canada gap-filling) 

	FracSEBClosure=1; fPlot=0; 

	[NEP,RE,GPP,NEPgf,REgf,GPPgf,REhat,GPPhat] = ...
		abNacpFcrnCO2Flux2NEP20090205 ...
		(t,NEE,uStar,Rsd,Ta,Ts,Rsd,Ta,Ts,fNight,uStarTh,FracSEBClosure,cSiteYr,fPlot); 
	
	NEPhat=GPPhat-REhat; xNEEhat=-NEPhat; xNEE=-NEP; 

%	=======================================================================

%	Characterize the random uncertainty curve 
% 	using the function abNacpRandomUncertaintyVPlot20100105, which: 
%	-	stratifies the year into NEEhat classes, 
%	-	computes the double expontential parameter mu for each NEEhat class, 
%	-	fits mu = f(NEEhat) lines for positive and negative values of NEEhat, 
%	-	estimates MuHat = f(NEEhat) for each period from the fitted lines.
	
	iFig=100; 
	
	[MuHat,b,cib,r2,s,mxNeg,myNeg,mxPos,myPos]= ...
		abNacpRandomUncertaintyVPlot20100105_2(xNEE,xNEEhat,iFig,cSiteYr,cFlux);
% 	b=[bNeg bPos];
    bNeg=b(:,1);
    bPos=b(:,2);
	if iFig>0; 
		eval(['print -djpeg100 ' DirOut 'randomUncertainty\figures\nacpAndrewPlot_' cSiteYr cFlux ';']);
	end; 
% 	eval(['save ' DirOut 'randomUncertainty\StatsNacpAndrewPlot_' cSiteYr ' b cib r2 s mxNeg myNeg mxPos myPos;']);
	
%	=======================================================================

% 	Random uncertainty analysis.

%	Initialize matrices for the aggregated Ensemble data products. 
%	-	Annual sum a*
%	-	Seasonal cycle at three integrating periods:
%		daily, weekly and monthly sums sd* sw* sm*.
%	-	Diurnal cycle at three averaging periods:
%		weekly, monthly and annual dw* dm* and da*.
		
	[y,m,d]=fcDateVec(t);
	nt=length(t); 
	nRecsPerDay=round(1/nanmedian(diff(t)));
    
     % this section was giving errors from using round, tjk August 2012
    if nt<=8785
        nRecsPerDay=24;
    else
        nRecsPerDay=48;
    end
    
	nd=nt/nRecsPerDay; nw=52; nm=12;
	nMinsPerDay=24*60/nRecsPerDay;
	%k=12*nMinsPerDay*60/1e6; % convert to g C m-2. 
	
	MatrixNEE=NaN*ones(nt,nMC);
	
	aMT=NaN*ones(nMC,1); aNEP=aMT; aRE=aMT; aGPP=aMT; % a annual
	nMissNEP=aMT; nMissRE=aMT; nMissGPP=aMT; % a annual
	
	dwMT=NaN*ones(nMC,nRecsPerDay,nw); dwNEP=dwMT; dwRE=dwMT; dwGPP=dwMT; % dw weekly mean diurnal cycle
	dmMT=NaN*ones(nMC,nRecsPerDay,nm); dmNEP=dmMT; dmRE=dmMT; dmGPP=dmMT; % dm monthly mean diurnal cycle
	daMT=NaN*ones(nMC,nRecsPerDay); daNEP=daMT; daRE=daMT; daGPP=daMT; % da annual mean diurnal cycle
	
	stMT=NaN*ones(nMC,nt); stNEP=stMT; stRE=stMT; stGPP=stMT; % sd hourly mean seasonal cycle
	sdMT=NaN*ones(nMC,nd); sdNEP=sdMT; sdRE=sdMT; sdGPP=sdMT; % sd daily mean seasonal cycle
	swMT=NaN*ones(nMC,nw); swNEP=swMT; swRE=swMT; swGPP=swMT; % sw weekly mean seasonal cycle
	smMT=NaN*ones(nMC,nm); smNEP=smMT; smRE=smMT; smGPP=smMT; % sm monthly mean seasonal cycle
	
	clear *MT;
	
	itw=1:(nRecsPerDay*nw*7); % indices for weekly averaging, includes 52*7 days only. 
	
%	Start with synthetic data and add gaps.
	
	xNEEgappy=xNEEhat; iEx=find(isnan(xNEE)); xNEEgappy(iEx)=NaN;
	
%	Implement Monte-Carlo process	with nMC independent draws. 
	
	for iMC=1:nMC;
		
		SqRt2=2^0.5;
		r=SqRt2*fcDoubleExpRnd(1,nt,1);
		xNEENoise=r.*MuHat;
		xNEEnoisy=xNEEgappy+xNEENoise;
		
		uStar=NaN*t; uStarTh=0; FracSEBClosure=1; % already dealt with
		
		[NEP,RE,GPP,NEPgf,REgf,GPPgf,REhat,GPPhat,REhat0,GPPhat0,mtRE,mcRE,mtGPP,mcGPP,bRE,bGPP] = ...
			abNacpFcrnCO2Flux2NEP20090205 ...
			(t,xNEEnoisy,uStar,Rsd,Ta,Ts,Rsd,Ta,Ts,fNight,uStarTh,FracSEBClosure,cSiteYr,fPlot); 
		
		%	Assign raw and aggregated outputs at different time scales and cycles.
		
		MatrixNEE(:,iMC)=-NEPgf;
		
		%	Integrated annual seasonal cycles, g C m-2 t-1.
		
		aNEP(iMC)=k*nansum(NEPgf); nMissNEP(iMC)=sum(isnan(NEPgf));
		aRE(iMC)=k*nansum(REgf); nMissRE(iMC)=sum(isnan(REgf));
		aGPP(iMC)=k*nansum(GPPgf); nMissGPP(iMC)=sum(isnan(GPPgf));
		
		stNEP(iMC,:)=k*NEPgf;
		stRE(iMC,:)=k*REgf;
		stGPP(iMC,:)=k*GPPgf;
	   
        sdNEP(iMC,:)=k*nansum(reshape(NEPgf,nRecsPerDay,nd));
		sdRE(iMC,:)=k*nansum(reshape(REgf,nRecsPerDay,nd));
		sdGPP(iMC,:)=k*nansum(reshape(GPPgf,nRecsPerDay,nd));
	
        swNEP(iMC,:)=k*nansum(reshape(NEPgf(itw),nRecsPerDay*7,nw));
		swRE(iMC,:)=k*nansum(reshape(REgf(itw),nRecsPerDay*7,nw));
		swGPP(iMC,:)=k*nansum(reshape(GPPgf(itw),nRecsPerDay*7,nw));
		
		for im=1:12;
			it=(find(m==im)); nit=length(it);
			smNEP(iMC,im)=k*nansum(NEPgf(it));
			smRE(iMC,im)=k*nansum(REgf(it));
			smGPP(iMC,im)=k*nansum(GPPgf(it));
		end;
		
		%	Mean diurnal cycles, umol m-2 s-1.
		
		daNEP(iMC,:,:)=nanmean(reshape(NEPgf,nRecsPerDay,nt/nRecsPerDay),2);
		daRE(iMC,:,:)=nanmean(reshape(REgf,nRecsPerDay,nt/nRecsPerDay),2);
		daGPP(iMC,:,:)=nanmean(reshape(GPPgf,nRecsPerDay,nt/nRecsPerDay),2);
		
		dwNEP(iMC,:,:)=nanmean(reshape(NEPgf(itw),nRecsPerDay,7,nw),2);
		dwRE(iMC,:,:)=nanmean(reshape(REgf(itw),nRecsPerDay,7,nw),2);
		dwGPP(iMC,:,:)=nanmean(reshape(GPPgf(itw),nRecsPerDay,7,nw),2);
		
		for im=1:12;
			it=(find(m==im)); nit=length(it);
			dmNEP(iMC,:,im)=nanmean(reshape(NEPgf(it),nRecsPerDay,nit/nRecsPerDay),2);
			dmRE(iMC,:,im)=nanmean(reshape(REgf(it),nRecsPerDay,nit/nRecsPerDay),2);
			dmGPP(iMC,:,im)=nanmean(reshape(GPPgf(it),nRecsPerDay,nit/nRecsPerDay),2);
		end;
		
		disp(sprintf('Monte Carlo random uncertainty %s  %g/%g  %5.3f   NEE RE GPP %4.0f %4.0f %4.0f   nMiss %4.0f %4.0f %4.0f', ...
			cSiteYr,iMC,nMC,uStarTh, -aNEP(iMC),aRE(iMC),aGPP(iMC), nMissNEP(iMC),nMissRE(iMC),nMissGPP(iMC)));
		
	end; % for iMC=1:nMC;
	
	Ensemble.NEE.nMissing = nMissNEP; Ensemble.RE.nMissing = nMissRE; Ensemble.GPP.nMissing = nMissGPP; 
	
	Ensemble.NEE.Annual = -aNEP; Ensemble.RE.Annual = aRE; Ensemble.GPP.Annual = aGPP; 
	Ensemble.NEE.Monthly = -smNEP; Ensemble.RE.Monthly = smRE; Ensemble.GPP.Monthly = smGPP; 
	Ensemble.NEE.Weekly = -swNEP; Ensemble.RE.Weekly = swRE; Ensemble.GPP.Weekly = swGPP; 
	Ensemble.NEE.Daily = -sdNEP; Ensemble.RE.Daily = sdRE; Ensemble.GPP.Daily = sdGPP; 
	Ensemble.NEE.Hourly = -stNEP; Ensemble.RE.Hourly = stRE; Ensemble.GPP.Hourly = stGPP; 
	
	Ensemble.NEE.DielByYr = -daNEP; Ensemble.RE.DielByYr = daRE; Ensemble.GPP.DielByYr =  daGPP; 
	Ensemble.NEE.DielByMo = -dmNEP; Ensemble.RE.DielByMo = dmRE; Ensemble.GPP.DielByMo = dmGPP; 
	Ensemble.NEE.DielByWeek = -dwNEP; Ensemble.RE.DielByWeek = dwRE; Ensemble.GPP.DielByWeek = dwGPP; 
	
%	=======================================================================
	
	if iFig>0 && strcmp(cFlux,'NEE'); 		
	
		fcFigLoc(iFig,0.5,0.45,'SE'); 
		
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
    elseif iFig>0
        % for LE and H don't plot components
        
        fcFigLoc(iFig,0.3,0.4,'NW'); 
		
		for iFlx=1:1; 
			
			switch iFlx; case 1; a=-aNEP; cFlx=[cSiteYr ' ' cFlux]; case 2; a=aRE; cFlx='RE'; case 3; a=aGPP; cFlx='GPP'; end; 
			subplot('position',[0.08+(iFlx-1)*0.32 0.76 0.76 0.18]); box on; 
			hist(a,30); title(sprintf('%s %g (%g)',cFlx,round(mean(a)),round(std(a)))); 
			if iFlx==1; ylabel('Histogram '); end; 
		
			switch iFlx; case 1; sm=-smNEP; cFlx=cFlux; case 2; sm=smRE; cFlx='RE'; case 3; sm=smGPP; cFlx='GPP'; end; 
			sp=[0.08+(iFlx-1)*0.32 0.41 0.76 0.28]; 
			subplot('position',sp); box on; 
			boxplot(sm); set(gca,'position',sp); 
			if iFlx==1; ylabel('Monthly '); end; 
			
			switch iFlx; case 1; da=-daNEP; case 2; da=daRE; case 3; da=daGPP; end; 
			sp=[0.08+(iFlx-1)*0.32 0.06 0.76 0.28]; 
			subplot('position',sp); box on; 
			boxplot(da); set(gca,'position',sp); 
			if iFlx==1; ylabel('Diel '); end; 
	
		end; 	
	end; 

%	=======================================================================
%	=======================================================================


