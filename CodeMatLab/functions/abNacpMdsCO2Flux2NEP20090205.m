function [xNEE,xNEEgf,xNEEhat] = ... 
	abNacpMdsCO2Flux2NEP20090205 ... 
	(t,NEE,uStar,PPFD,Ta,Ts,Vpd,PPFDGF,TaGF,TsGF,VpdGF,isNight,uStarTh,iFig,cSiteYr); 

%abNacpMdsCO2Flux2NEP20090205
%
%FC created 29 Sept 2010 to reset to normal FCRN program. 
%but with small filled gaps assigned to gf not measured fluxes. 
%
%makes some slight modifications to FCRN_CO2Flux2NEP
%-uses nlinfit which is faster
%-deals with large gaps
%-stripped down for speed. 
%-uses ratio not forced-origin lin reg to estimate mcR and mcGEP
%-reduces Inc to 10.

%-	post-processes eddy-covariance measurements of 
%	net ecosystem exchange (NEE = CO2 flux + storage)
%-	calculates net ecosystem productivity (NEP) 
%-	partitions NEP into ecosystem respiration (R) 
%	and gross ecosystem photosynthesis (GEP)
%-	fills gaps in NEP, R and GEP, producing NEPgf, Rgf, and GEPgf
%
%Processing Steps:
%----------------
%
%1. NEP is estimated as -NEE after excluding low u* data at night and
%	applying the optional energy closure adjustment.
%2. Measured R is estimated from NEE at night and during the non-GrowingSeason
%3. An empirical R=f(Ts,t) model is fit based on the data from 2
%4. The R=f(Ts,t) model from 3 is used to estimate R during the day 
%	and to fill gaps in R at night
%5. GEP is estimated as NEP+R (daytime, growing season) or zero 
%	(nighttime and non-GrowingSeason) 
%6. An empirical GEP=f(PPFD,t) model is fit to the GEP data from 5
%7. The model from 6 is used to fill gaps in daytime, growing-season GEP
%8. Gaps in NEP are filled using modelled GEP-R
%
%	A moving window is used to capture the time variation in 
%	the R and GEP model parameters cR(t) and cGEP(t): 
%		R 		= cR(t)*bR(1)/(1+exp[bR(2)*(bR(3)-Ts)]); 
%		GEP 	= cGEP(t)*bGEP(1)*PPFD/(bGEP(2)+PPFD); 
%
%Syntax:
%------
%
%[NEP,R,GEP,NEPgf,Rgf,GEPgf] 
%	= CO2Flux2NEP(t,NEE,uStar,PPFD,Ta,Ts,PPFDGF,TaGF,TsGF,isNight,uStarTh,FracSEBClosure); 
%
%Input Arguments:
%---------------
%
%-	t is the decimal day vector
%-	NEE is the measured NEE (CO2 flux + storage) vector
%-	uStar is the friction velocity vector
%-	PPFD is the downwelling PAR flux density vector (used to model GEP)
%-	Ta is the air temperature vector (used optionally to model GEP)
%-	Ts is the soil temperature vector (used to model R)
%-	PPFDGF, TaGF, and TsGF are gap-filled values of PPFD, Ta, and Ts  
%-	isNight is a vector that indicates whether itReg is day (0) or night (1)
%-	uStarTh is the scalar (or vector) uStar threshold below which nighttime NEP data are rejected
%-	FracSEBClosure is a scalar or vector SEB-closure adjustment divider 
%-	LogFile is the (optional) log file directory and name; create if not empty. 
%-	Plots (optional) specifies the plots to create. 

%	========================================================================
%	========================================================================

%	Notes

%	FCRN_CO2Flux2NEP Version 1.0
%	Written by Alan Barr and Natascha Kljun, 28 July 2003
%	Last update, 12 Sept 2003 FROZEN AS FCRN VERSION 1.0. 

%	Several non-standard MatLab functions are celled: 
%	- doy
%	- FillSmallGapsByLinInterp 
%	- fr_function_fit
%	- gapsize 
%	- mydatevec
%	- mynanmean
%	- mynanmedian
%	- mynansum
%	- naninterp1

%	In this implementation, the growing season is not preset externally but falls 
%	out of the GEP analysis (in the cGEP(t) parameter). It is also constrained by air
%	and soil temperatures, which are used to stratify that data into three seasons:
%	-	the cold season (when both Ta and Ts are less than or equal to zero), 
%		when daytime GEP is set to zero, 
%	-	the non-cold season (when either or both Ta and Ts are positive), 
%		when daytime GEP is estimated from measured NEP and modeled R, 
%	-	the warm season (when both Ta and Ts are positive), which is used to estimate 
%		the bGEP model parameters.

%	========================================================================
%	========================================================================

	DirCurrent=pwd; 
	DirMds='i:\nacpUncertaintyExamples\CodeMDS\'; 
	cd(DirMds); 
	
	cSite=cSiteYr(1:5); iYr=str2num(cSiteYr(7:10)); 

%	========================================================================

%	Exclude low-u* data at night. 

	iNight = find(isNight); 
	iLowuStar = find(uStar<uStarTh); 
	iExclude = intersect(iLowuStar,iNight); 
	NEE(iExclude) = NaN; 
	
%	Implement MDS gap filling in three steps: 
%	-	write ascii input file
%	-	call gapfilling executable 
%	-	read and assign gap-filling output

	Dir='';
	NamesMds={'NEE','RG','TA','VPD'};
	xMds=[NEE PPFDGF TaGF VpdGF];
	
	FileMdsIp=nacpWriteMdsInputFile(cSite,Dir,t,xMds,NamesMds); % missing -9999
	
	eval(['!gapfilling -i' FileMdsIp ' -e']); % The exclamation mark ! invokes an operating system comand.
	
	[xGF,NamesGF,fErr]=nacpReadMdsOutputFile(DirMds,cSite,iYr,1);
	
	if fErr==1;
		[xGF,NamesGF,fErr]=nacpReadOldMdsOutputFile(DirMds,cSite,iYr,1);
	end;
	if fErr;
		nacpAddWhitespaceToMdsOutputFile(DirMds,cSite,iYr,1);
		[xGF,NamesGF]=nacpReadMdsOutputFile(DirMds,cSite,iYr,1);
	end;
	FileMdsOp=strrep(FileMdsIp,'.txt','_gap.txt');
	
	icNEE=strmatch('ORIG_NEE',NamesGF,'exact'); xNEE=xGF(:,icNEE);
	icNEEhat=strmatch('FILLED_NEE',NamesGF,'exact'); xNEEhat=xGF(:,icNEEhat);
	
	xNEEgf=xNEE; iGF=find(isnan(xNEEgf) & ~isnan(xNEEhat)); xNEEgf(iGF)=xNEEhat(iGF);
	
	cd(DirCurrent); 
	
%	========================================================================
%	========================================================================

	if iFig>0; 
		
		myFigLoc(iFig,0.3,0.35,'SE');
		
		plot(t,xNEEgf,'r.',t,xNEE,'b.'); mydatetick(t,'Mo',4,1); xlim([min(t) max(t)]);
		title(sprintf('%s  NEE nMiss %g',cSiteYr,sum(isnan(myrv(xNEEgf)))));
		box on; 
		
	end;
	
%	========================================================================
%	========================================================================

	