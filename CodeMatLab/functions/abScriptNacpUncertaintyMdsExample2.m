%	abScriptNacpUncertaintyMdsExample2

%	is a sample script that provides an example of Alan Barr's 
%	2010 implementation of the NACP NEE uncertainy analysis 
%	for some typical FLUXNET L3 files, 
%	using Markus Reichstein's MDS gap filling.  

%	========================================================================
%	========================================================================

%	In-house functions called: 

%	abCpdBootstrapUStarTh20100901
%	abCpdAssignUStarTh20100901
%	abNacpMdsCO2Flux2NEP20090205
%	abNacpMdsUStarTh2NepUncertainty20100423
%	abNacpMdsRandomUncertainty20100423
%	abNacpNeeOutlierExclusion
%	abNacpPlotUncertaintyInputs
%	fcx2rowvec
%	fcEqnAnnualSine
%	fcReadFields
%	fcr2Calc
%	fcFigLoc

%	stats toolbox:  boxplot, nanmedian, nanmean, nlinfit

%	Written 16 April 2010 by Alan Barr

%	========================================================================
%	========================================================================

	close all; clear all; clc; warning off; 
	
	addpath('i:\nacpUncertaintyExamples\CodeMatLab\'); 
		
	DirRoot='i:\nacpUncertaintyExamples\'; 
	
	mkdir([DirRoot 'Output\MDS\']); DirOut=[DirRoot 'Output\MDS\']; 
	mkdir([DirRoot 'Data_L3\']); DirL3=[DirRoot 'Data_L3\']; 
	mkdir([DirRoot 'figuresQC\']); DirQC=[DirRoot 'figuresQC\']; 
	mkdir([DirRoot 'log\']); DirLog=[DirRoot 'Log\']; 

	mkdir([DirOut 'randomUncertainty\']); 
	mkdir([DirOut 'uStarThAnalysis\']); 
	mkdir([DirOut 'uStarTh2NeeUncertainty\']); 
	
	mkdir([DirOut 'randomUncertainty\figures\']); 
	mkdir([DirOut 'uStarThAnalysis\figures\']); 
	mkdir([DirOut 'uStarTh2NeeUncertainty\figures\']); 

%	========================================================================

%	Input one year of data at a time. 

	Files=dir([DirL3 '*_L3.mat']); nFiles=length(Files); 
	
	for iFile=1:nFiles; 
		
		FileL3=Files(iFile).name; 
		i=strfind(FileL3,'_L3');
		ii=(i-9):(i-5); cSite=FileL3(ii); 
		ii=(i-4):(i-1); cYr=FileL3(ii); iYr=str2num(cYr); 
		cSiteYr=[cSite '-' cYr]; 
		
		disp(sprintf('abScriptNacpUncertaintyMdsExample2 is processing file %g of %g, %s',iFile,nFiles,cSiteYr)); 
		disp(datestr(now)); 
		disp(' ' ); 
	
		FileDiary=[DirLog cSiteYr '.log']; diary(FileDiary); 

%	load

	load ([DirL3 FileL3]); 
	
	data_L3(data_L3==-9999)=NaN; % set missing to NaN
	data_L3(data_L3==-6999)=NaN; % set missing to NaN
	
	i=strmatch('DoY',int_L3,'exact'); d=data_L3(:,i); 
	i=strmatch('ustar',int_L3,'exact'); uStar=data_L3(:,i); 
	i=strmatch('NEE_or',int_L3,'exact'); NEE=data_L3(:,i); 
	i=strmatch('Ta',int_L3,'exact'); Ta=data_L3(:,i); 
	i=strmatch('Ts1',int_L3,'exact'); Ts=data_L3(:,i); 
	i=strmatch('PPFD',int_L3,'exact'); PPFD=data_L3(:,i); 
	i=strmatch('Rg',int_L3,'exact'); Rg=data_L3(:,i); 
	i=strmatch('Rh',int_L3,'exact'); Rh=data_L3(:,i);

	if sum(~isnan(NEE))==0; 
		i=strmatch('NEE_st',int_L3,'exact'); NEE=data_L3(:,i);
	end; 	
	
	t=datenum(iYr-1,12,31)+d; nt=length(t); 
	nRecsPerDay=round(1/nanmedian(diff(t)));
	
	es=10*fcekPaTetenAboveWater(Ta,100); 
	Vpd=es.*(1-Rh/100); Vpd(find(Vpd<0))=0; 
	VpdGF=Vpd; 

	fNight=Rg<5; % flag nighttime periods
	fNight=PPFD<10; % flag nighttime periods

%	========================================================================
	
%	Mildly filter NEE to exclude extreme outliers

	nsP=7; nsCi=7; ndWindow=28; iFig=201; 

	[NeeQC,iOut] = abNacpNeeOutlierExclusion ... 
		(t,NEE,PPFD,Ta,Ts,fNight,nsP,nsCi,ndWindow,cSiteYr,iFig); 
	nOut=length(iOut);
	
	disp(' ');
	disp(sprintf([cSiteYr ' %4.0f outliers excluded.'], nOut));
	disp(' ');

	if iFig>0; 
		eval(['print -djpeg100 ' DirQC 'NeeOutliers_' cSiteYr ';']); 
	end; 
	
	NEE=NeeQC; 
	
%	========================================================================
	
%	Plot input data for visual inspection. 

	iFig=202; 
	
	abNacpPlotUncertaintyInputs(t,NEE,uStar,PPFD,Ta,Ts,Vpd,fNight,cSiteYr,iFig); 
	
	if iFig>0; 
		eval(['print -djpeg100 ' DirQC 'PlotInputs_' cSiteYr ';']); 
	end;  
			
%	========================================================================

%	uStarTh analysis

%	Call uStarTh bootstrappng program and assign annual Cp (change-point) arrays.	
	
	nBoot=1e2; iFig=203; fPlot=0; 
	
	[Cp2,Stats2,Cp3,Stats3] = ... 
		abCpdBootstrapUStarTh20100901 ...
			(t,NEE,uStar,Ts,fNight,fPlot,cSiteYr,nBoot); 
	[CpBoot,n,tW,CpW,cMode,cFailure,fSelect,sSine,FracSig,FracModeD,FracSelect,CpArray] ... 
		= abCpdAssignUStarTh20100901(Stats2,iFig,cSiteYr); 
	
%	Assign various Cp (the change-point, alternate designation for uStarTh) 	
	
	Cp=mean(CpBoot); 
	CpBoot=fcx2rowvec(CpBoot); 
	CpSine=fcEqnAnnualSine(sSine,t); 
	
	%	assign seasonal uStarTh sine curve for each bootstrap (CpSineBoot); 
	
	mt=fcReadFields(Stats2,'mt'); [nW,nS,nB]=size(mt); 
	mt=reshape(mt,nW*nS,nB); CpArray=reshape(CpArray,nW*nS,nB); 
	CpSineBoot=NaN*ones(nt,nBoot); 
	for i=1:nBoot; 
		bSine=[1,1,1]; xx=mt(:,i); yy=CpArray(:,i); 
		iNaN=find(isnan(xx+yy)); xx(iNaN)=[]; yy(iNaN)=[]; 
		[bSine]=nlinfit(xx,yy,'fcEqnAnnualSine',bSine); 
		yHat=fcEqnAnnualSine(bSine,xx); r2=fcr2Calc(yy,yHat); 
		if bSine(2)<0; bSine(2)=-bSine(2); bSine(3)=bSine(3)+365.2425/2; end; 
		bSine(3)=mod(bSine(3),365.2425); 
		b1=bSine(1); b2=bSine(2); b3=bSine(3); 
		CpSineBoot(:,i)=fcEqnAnnualSine([b1 b2 b3 r2],t); 
	end; 
	
	eval(['save ' DirOut 'uStarThAnalysis\uStarThAnalysis_' cSiteYr ... 
		' cSite cYr CpBoot CpSine CpSineBoot n tW CpW cMode cFailure fSelect sSine FracSig FracModeD FracSelect;']);

	if iFig>0; 
		eval(['print -djpeg100 ' DirOut 'uStarThAnalysis\figures\uStarTh_' cSiteYr ';']); 
	end; 
	
%	Decide whether to use single annual or seasonal sine uStarTh. 
%	Could be based on goodness of fit of anual sine. 

	if sSine(2)./sSine(1)>0.2 & sSine(4)>0.2;
			whichCp='Seasonal';
		else;
			whichCp='Annual';
	end;
	
	disp(' '); disp(['Using ' whichCp ' uStarTh.']); disp(' '); 
	
	switch whichCp; 
		case 'Annual'; uStarThs=CpBoot;
		case 'Seasonal'; uStarThs=CpSineBoot; 
			fcFigLoc(99,0.30,0.35,'NC'); 
			nm=round(nt/12); imt=(nm/2):nm:nt; 
			boxplot(uStarThs(imt,:)'); box on; 
			xlabel('Month'); ylabel('u_*^{Th} (m s^{-1})'); 
			title(sprintf('%s   s_0=%4.2f  s_1=%4.2f  r^2=%4.2f ',cSiteYr,sSine([1 2 4]))); 
			eval(['print -djpeg100 ' DirOut 'uStarThAnalysis\figures\uStarThsSineBoot_' cSiteYr ';']); 
	end; 
		
%	========================================================================

%	check processing at mean uStarTh. 

	FracSEBClosure=1; muStarTh=mean(uStarThs,2); iFig=204; 
		
	[xNEE,xNEEgf,xNEEhat] = ...
		abNacpMdsCO2Flux2NEP20090205 ...
		(t,NEE,uStar,PPFD,Ta,Ts,Vpd,PPFD,Ta,Ts,Vpd,fNight,muStarTh,iFig,cSiteYr);
	
	NEPgf=-xNEEgf; REgf=NaN*t; GPPgf=NaN*t; 
	
	nMinsPerDay=24*60/nRecsPerDay; k=12*nMinsPerDay*60/1e6; 
	
	aNEP=k*nansum(NEPgf); nNEP=sum(~isnan(NEPgf)); nMissNEP=sum(isnan(NEPgf));
	aRE=k*nansum(REgf); nMissRE=sum(isnan(REgf));
	aGPP=k*nansum(GPPgf); nMissGPP=sum(isnan(GPPgf));
	
	disp(sprintf('Filling gaps in %s, u*Th %5.3f   %4.0f %4.0f %4.0f   n %5.0f nMiss %4.0f %4.0f %4.0f', ...
		cSiteYr, mean(muStarTh), aNEP,aRE,aGPP, nNEP, nMissNEP,nMissRE,nMissGPP)); disp(' '); 
	
	if iFig>0; 
		eval(['print -djpeg100 ' DirQC 'GapFillingMds_' cSiteYr ';']); 
	end;

%	========================================================================

% 	calc uStarTh-related uncertainties in NEE, RE and GPP

	iFig=205; 
	
	[uuEnsemble,uuNEE] = ... 
		abNacpMdsUStarTh2NepUncertainty20100423 ... 
		(t,NEE,uStar,PPFD,Ta,Ts,Vpd,PPFD,Ta,Ts,Vpd,fNight,uStarThs,iFig,cSiteYr); 
	
	eval(['save ' DirOut 'uStarTh2NeeUncertainty\uStarTh2NeeMdsUncertainty_' cSiteYr ... 
		' cSite cYr uuEnsemble,uuNEE;']);
	
	if iFig>0; 
		eval(['print -djpeg100 ' DirOut 'uStarTh2NeeUncertainty\figures\uStarTh2NeeMdsUncertainty_' cSiteYr ';']); 
	end; 
	
%	========================================================================

%	Random uncertainty analysis

%	Package as function separately for NEE and GPP RE
%	have entire array output as an option. 

% % % 	uStarTh=mean(CpBoot); 

	uStarTh=CpSine; nMC=1e2; iFig=206; % nMC is the # MonteCarlo reps

	[ruEnsemble,ruNEE]=abNacpMdsRandomUncertainty20100423 ...
		(t,NEE,uStar,PPFD,Ta,Ts,Vpd,fNight,uStarTh,DirOut,cSiteYr,nMC,iFig); 
		
	eval(['save ' DirOut 'randomUncertainty\randomUncertaintyMds_' cSiteYr ... 
		' cSite cYr ruEnsemble ruNEE;']);
 
	if iFig>0; 
		eval(['print -djpeg100 ' DirOut 'randomUncertainty\figures\randomUncertaintyMds_' cSiteYr ';']); 
	end; 
	
%	========================================================================

	disp(' '); disp(' '); 

	diary off; 
	
	end; 
	
	