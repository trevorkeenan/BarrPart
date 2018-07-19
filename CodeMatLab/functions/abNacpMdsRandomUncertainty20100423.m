	function [Ensemble,MatrixNEE]=abNacpMdsRandomUncertainty20100423 ...
		(t,NEE,uStar,Rsd,Ta,Ts,Vpd,fNight,uStarTh,DirOut,cSiteYr,nMC,iFig) 


%abNacpMdsRandomUncertainty20100423
%	random uncertainty analysis function, 
%	based on Richardson & Hollinger (2007), 
%	using Markus Reichstein's MDS
%	gap-filling procedure,  
%	typically run one year at a time. 
%
%Syntax: 
%
%	[Ensemble,MatrixNEE] = abNacpMdsRandomUncertainty20100423 ...
% 		(t,NEE,uStar,Rsd,Ta,Ts,Vpd,fNight,uStarTh,DirOut,cSiteYr,nMC); 
%
%	-	Ensemble is a structured record containing the ensemble output
%		(typically nMC = 1,000 realizations) of NEE, 
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

	[xNEE,xNEEgf,xNEEhat] = ...
		abNacpMdsCO2Flux2NEP20090205 ...
		(t,NEE,uStar,Rsd,Ta,Ts,Vpd,Rsd,Ta,Ts,Vpd,fNight,uStarTh,fPlot,cSiteYr); 
	
%	=======================================================================

%	Characterize the random uncertainty curve 
% 	using the function abNacpRandomUncertaintyVPlot20100105, which: 
%	-	stratifies the year into NEEhat classes, 
%	-	computes the double expontential parameter mu for each NEEhat class, 
%	-	fits mu = f(NEEhat) lines for positive and negative values of NEEhat, 
%	-	estimates MuHat = f(NEEhat) for each period from the fitted lines.
	
	iFig=100; 
	
	[MuHat,b,cib,r2,s,mxNeg,myNeg,mxPos,myPos]= ...
		abNacpRandomUncertaintyVPlot20100105(xNEE,xNEEhat,iFig,cSiteYr);
	
	if iFig>0; 
		eval(['print -djpeg100 ' DirOut 'randomUncertainty\figures\nacpAndrewPlot_' cSiteYr ';']);
	end; 
	eval(['save ' DirOut 'randomUncertainty\StatsNacpAndrewPlot_' cSiteYr ' b cib r2 s mxNeg myNeg mxPos myPos;']);
	
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
	nd=nt/nRecsPerDay; nw=52; nm=12;
	nMinsPerDay=24*60/nRecsPerDay;
	k=12*nMinsPerDay*60/1e6; % convert to g C m-2. 
	
	MatrixNEE=NaN*ones(nt,nMC);
	
	aMT=NaN*ones(nMC,1); aNEP=aMT; % a annual
	nMissNEP=aMT; 
	
	dwMT=NaN*ones(nMC,nRecsPerDay,nw); dwNEP=dwMT; % dw weekly mean diurnal cycle
	dmMT=NaN*ones(nMC,nRecsPerDay,nm); dmNEP=dmMT; % dm monthly mean diurnal cycle
	daMT=NaN*ones(nMC,nRecsPerDay); daNEP=daMT; % da annual mean diurnal cycle
	
	sdMT=NaN*ones(nMC,nd); sdNEP=sdMT; % sd daily mean seasonal cycle
	swMT=NaN*ones(nMC,nw); swNEP=swMT; % sw weekly mean seasonal cycle
	smMT=NaN*ones(nMC,nm); smNEP=smMT; % sm monthly mean seasonal cycle
	
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
		
		[mcNEE,mcNEEgf] = abNacpMdsCO2Flux2NEP20090205 ...
			(t,xNEEnoisy,uStar,Rsd,Ta,Ts,Vpd,Rsd,Ta,Ts,Vpd,fNight,uStarTh,fPlot,cSiteYr); 
		
		%	Assign raw and aggregated outputs at different time scales and cycles.
		
		MatrixNEE(:,iMC)=mcNEEgf; NEPgf=-mcNEEgf; 
		
		%	Integrated annual seasonal cycles, g C m-2 t-1.
		
		aNEP(iMC)=k*nansum(NEPgf); nMissNEP(iMC)=sum(isnan(NEPgf));
		sdNEP(iMC,:)=k*nansum(reshape(NEPgf,nRecsPerDay,nd));
		swNEP(iMC,:)=k*nansum(reshape(NEPgf(itw),nRecsPerDay*7,nw));
		for im=1:12;
			it=(find(m==im)); nit=length(it);
			smNEP(iMC,im)=k*nansum(NEPgf(it));
		end;
		
		%	Mean diurnal cycles, umol m-2 s-1.
		
		daNEP(iMC,:,:)=nanmean(reshape(NEPgf,nRecsPerDay,nt/nRecsPerDay),2);
		dwNEP(iMC,:,:)=nanmean(reshape(NEPgf(itw),nRecsPerDay,7,nw),2);
		for im=1:12;
			it=(find(m==im)); nit=length(it);
			dmNEP(iMC,:,im)=nanmean(reshape(NEPgf(it),nRecsPerDay,nit/nRecsPerDay),2);
		end;
		
		disp(sprintf('Monte Carlo random uncertainty %s  %g/%g  %5.3f   NEE %4.0f   nMiss %4.0f ', ...
			cSiteYr,iMC,nMC,uStarTh, -aNEP(iMC), nMissNEP(iMC)));
		
	end; % for iMC=1:nMC;
	
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


