function [ruEnsemble, uuEnsemble, NeeQC,NEEgfgC,bNeg,bPos] ...
    = abScriptNacpUncertaintyFcrn_fixedUstarDP(dataL2_4alan,cSite,cYr,columns,nBoot,uThresh)

%	abScriptNacpUncertaintyFcrnExample1

%	is a sample script that provides an example of Alan Barr's
%	2010 implementation of the NACP NEE uncertainy analysis
%	for a sample L2 file from Trevor Keenan.

%	program compromises to run this year
%	-	Fc used instead of NEE, which was u*-filtered
%	-	relaxed data requirements for uStarTh
%	-	uStarTh-related uncertainty based on annual sine fit
%		to avoid excessive exclusion of summer data.

%	========================================================================
%	========================================================================

%	In-house functions called:

%	abCpdBootstrapUStarTh20100901
%	abCpdAssignUStarTh20100901
%	abNacpFcrnCO2Flux2NEP20090205
%	abNacpFcrnUStarTh2NepUncertainty20100423
%	abNacpFcrnRandomUncertainty20100423
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

cFlux='NEE';
   
close all;
%     clear all;

warning off;

% addpath('./functions/');
% addpath('./functions2/')


mkdir('../../Folder-Structure_Results/Figures/'); 
DirRoot='../../Folder-Structure_Results/Figures/';
DirRoot2='../../Folder-Structure_Results/Figures/';

mkdir([DirRoot2 'Output/FCRN/']); DirOut=[DirRoot2 'Output/FCRN/'];
% mkdir([DirRoot 'Data_L2/']); %DirL2=[DirRoot 'Data_L2/'];
mkdir([DirRoot2 'figuresQC/']); DirQC=[DirRoot2 'figuresQC/'];
mkdir([DirRoot2 'log/']); DirLog=[DirRoot2 'Log/'];


%	========================================================================

cSiteYr=[cSite '-' cYr];

FileDiary=[DirLog cSiteYr '.log']; diary(FileDiary);


% 	assign
year=dataL2_4alan(:,columns.year);
ddoy=dataL2_4alan(:,columns.ddoy);

% 	Fc=dataL2_4alan(:,columns.Fc);
NEE=dataL2_4alan(:,columns.NEE);

uStar=dataL2_4alan(:,columns.ustar);
Ta=dataL2_4alan(:,columns.Tair);
Ts=dataL2_4alan(:,columns.SoilT);
PPFD=dataL2_4alan(:,columns.PAR);
Vpd=dataL2_4alan(:,columns.VPD);

NEE(NEE<=-999)=NaN;
uStar(uStar==-9999)=NaN;
Ta(Ta==-9999)=NaN;
Ts(Ts==-9999)=NaN;
PPFD(PPFD==-9999)=NaN;
Vpd(Vpd==-9999)=NaN;


% I do not have Ts
% Set = Ta
if sum(~isnan(Ts))==0
    Ts=Ta;
end
% as we don't have Tsoil
% set soil T to zero at start and end of year to mark a cold season
% Ts(1:24*30)=0;
% Ts(end-24*30:end)=0;


%	derive and accomodate
t=datenum(year,ones(length(year),1),ddoy); nt=length(t);
nRecsPerDay=round(1/nanmedian(diff(t)));
% if nt >1700
%     nRecsPerDay=48;
% else
%     nRecsPerDay=24;
% end
nMinsPerDay=24*60/nRecsPerDay;
k=12*nMinsPerDay*60/1e6;

% change here to change the night time detection limit
% Ts should be Tsoil 0 need to use gap filled Tsoil, tjk
fNight=PPFD<20; % had to relax from 10 to 20
% 	Ts=Ta; % Ts has gaps; prefer to use shallow Ts or mean of Ta and Ts.

% 	NEE=Fc;	% bogus temporary fix because raw NEE is not included
%	The included NEE is u*-fitered and gap-filled.

%	========================================================================

%	Mildly filter NEE to exclude extreme outliers

% 	nsP=7; nsCi=7; ndWindow=28; iFig=101;
nsP=7; nsCi=7; ndWindow=20; iFig=0;

[NeeQC,iOut] = abNacpNeeOutlierExclusion_2 ...
    (t,NEE,PPFD,Ta,Ts,fNight,nsP,nsCi,ndWindow,cSiteYr,iFig,cFlux);
nOut=length(iOut);

disp(' ');
disp(sprintf([cSiteYr ' %4.0f outliers excluded.'], nOut));
disp(' ');

if iFig>0;
    eval(['print -djpeg100 ' DirQC 'NeeOutliers_' cSiteYr ';']);
end;

NEE=NeeQC;
NeeQC=k*NeeQC;

%	========================================================================

%	Plot input data for visual inspection.

iFig=102;

abNacpPlotUncertaintyInputs(t,NEE,uStar,PPFD,Ta,Ts,Vpd,fNight,cSiteYr,iFig);

if iFig>0;
    eval(['print -djpeg100 ' DirQC 'PlotInputs_' cSiteYr ';']);
end;

%	========================================================================

%	uStarTh analysis

%	Call uStarTh bootstrappng program and assign annual Cp (change-point) arrays.

iFig=103; fPlot=0;

% nboot usually set to 1000 but will take an hour to run

% the following function runs the bootstrapping, tjk
% [Cp2,Stats2,Cp3,Stats3] = ...
%     abCpdBootstrapUStarTh20100901 ...
%     (t,NEE,uStar,Ts,fNight,fPlot,cSiteYr,nBoot);
% 
% % this function does quality assurance of the ustar change point,
% % tjk
% [CpBoot,n,tW,CpW,cMode,cFailure,fSelect,sSine,FracSig,FracModeD,FracSelect,CpArray] ...
%     = abCpdAssignUStarTh20100901(Stats2,iFig,cSiteYr);

%	Assign various Cp (the change-point, alternate designation for uStarTh)
% print cFailure to screen

% cFailure
% Cp=nanmean(CpBoot);

% tjk, Aug 2012
CpBoot=uThresh;
Cp=uThresh;

% if ~isnan(Cp)   % if we don't find a Cp, then skip all
%     CpBoot=fcx2rowvec(CpBoot);
%      CpSine=fcEqnAnnualSine(sSine,t);
%     
%     %	assign seasonal uStarTh sine curve for each bootstrap (CpSineBoot);
%     
%     mt=fcReadFields(Stats2,'mt'); [nW,nS,nB]=size(mt);
%     mt=reshape(mt,nW*nS,nB); CpArray=reshape(CpArray,nW*nS,nB);
%     CpSineBoot=NaN*ones(nt,nBoot);
%     for i=1:nBoot;
%         bSine=[1,1,1]; xx=mt(:,i); yy=CpArray(:,i);
%         iNaN=find(isnan(xx+yy)); xx(iNaN)=[]; yy(iNaN)=[];
%         [bSine]=nlinfit(xx,yy,'fcEqnAnnualSine',bSine);
%         yHat=fcEqnAnnualSine(bSine,xx); r2=fcr2Calc(yy,yHat);
%         if bSine(2)<0; bSine(2)=-bSine(2); bSine(3)=bSine(3)+365.2425/2; end;
%         bSine(3)=mod(bSine(3),365.2425);
%         b1=bSine(1); b2=bSine(2); b3=bSine(3);
%         CpSineBoot(:,i)=fcEqnAnnualSine([b1 b2 b3 r2],t);
%     end;
%     
%     eval(['save ' DirOut 'uStarThAnalysis/uStarThAnalysis_' cSiteYr ...
%         ' cSite cYr CpBoot CpSine CpSineBoot n tW CpW cMode cFailure fSelect sSine FracSig FracModeD FracSelect;']);
%     
%     if iFig>0;
%         eval(['print -djpeg100 ' DirOut 'uStarThAnalysis/figures/uStarTh_' cSiteYr ';']);
%     end;
%     
%     %	Decide whether to use single annual or seasonal sine uStarTh.
%     %	Could be based on goodness of fit of anual sine.
%     
%     if sSine(2)/sSine(1)>0.2 && sSine(4)>0.2;
%         whichCp='Seasonal';
%     else
%         whichCp='Annual';
%     end
%     
%     disp(' '); disp(['Using ' whichCp ' uStarTh.']); disp(' ');
%     
%     switch whichCp;
%         case 'Annual'; uStarThs=CpBoot;
%         case 'Seasonal'; uStarThs=CpSineBoot;
%             fcFigLoc(99,0.30,0.35,'NC');
%             nm=round(nt/12); imt=(nm/2):nm:nt;
%             boxplot(uStarThs(imt,:)'); box on;
%             xlabel('Month'); ylabel('u_*^{Th} (m s^{-1})');
%             title(sprintf('%s   s_0=%4.2f  s_1=%4.2f  r^2=%4.2f ',cSiteYr,sSine([1 2 4])));
%             eval(['print -djpeg100 ' DirOut 'uStarThAnalysis/figures/uStarThsSineBoot_' cSiteYr ';']);
%     end;
    
    %	========================================================================
    
    %	check processing at mean uStarTh.
    
    
uStarThs = Cp;  % tjk, Aug 2012

if nansum(NEE)~=0 && ~isnan(nansum(NEE))  % to deal with Howland not having any 2009 data
    FracSEBClosure=1; muStarTh=mean(uStarThs,2); iFig=110;%iFig=104;
        
    % gap fill the data once and plot the gap-filling figure at mean ustar
    % threshold, tjk
    [NEP,RE,GPP,NEPgf,REgf,GPPgf,RHat,GPPHat,RHat0,GPPHat0,mtR,mcR,mtGPP,mcGPP,bR,bGPP] = ...
        abNacpFcrnCO2Flux2NEP20090205 ...
        (t,NEE,uStar,PPFD,Ta,Ts,PPFD,Ta,Ts,fNight,muStarTh,FracSEBClosure,cSiteYr,iFig);
    
    
    aNEP=k*nansum(NEPgf); nNEP=sum(~isnan(NEPgf)); nMissNEP=sum(isnan(NEPgf));
    aRE=k*nansum(REgf); nMissRE=sum(isnan(REgf));
    aGPP=k*nansum(GPPgf); nMissGPP=sum(isnan(GPPgf));
    
    NEEgfgC=k*NEPgf;
    
    disp(sprintf('Filling gaps in %s, u*Th %5.3f   %4.0f %4.0f %4.0f   n %5.0f nMiss %4.0f %4.0f %4.0f', ...
        cSiteYr, mean(muStarTh), aNEP,aRE,aGPP, nNEP, nMissNEP,nMissRE,nMissGPP)); disp(' ');
    
    if iFig>0;
        eval(['print -djpeg100 ' DirQC 'GapFillingFcrn_' cSiteYr ';']);
    end;
    
    %	========================================================================
    
    % 	calc uStarTh-related uncertainties in NEE, RE and GPP
    
    iFig=0;
    
    % this outputs a structured record - uuEnsemble.NEE
    % where ustar uncertainty is stored
    % and ruEnsemble.NEE where random uncertainty is stored
    
    % check nMissing - should be full of zeros
    % if there's a big gap that will be highlighted here
    % maximum gap is 3*31 days.
    
    % take daily uncertainty with a grain of salt.
    % designed to work better at longer time scales
    
    % uuNEE is the whole NEE output for each bootstrap
    % RE and GPP could be output or aggregated also
    % just find the function where NEE is being done
    % and put in the RE and GPP stuff
    
    % stored in uStarTh2NeeFcrnUncertainty_UNkwn-9999
    
    
    
    [uuEnsemble,uuNEE] = ...
        abNacpFcrnUStarTh2NepUncertainty20100423 ...
        (t,NEE,uStar,PPFD,Ta,Ts,PPFD,Ta,Ts,fNight,uStarThs,cSiteYr,iFig,k);
    
%     eval(['save ' DirOut 'uStarTh2NeeUncertainty/uStarTh2NeeFcrnUncertainty_' cSiteYr ...
%         ' cSite cYr uuEnsemble,uuNEE;']);
     
    if iFig>0;
        eval(['print -djpeg100 ' DirOut 'uStarTh2NeeUncertainty/figures/uStarTh2NeeFcrnUncertainty_' cSiteYr ';']);
    end
    
    
    
    %	========================================================================
    
    %	Random uncertainty analysis
    
    %	Package as function separately for NEE and GPP RE
    %	have entire array output as an option.
    
    % % % 	uStarTh=mean(CpBoot);
    
    
    %uStarTh=CpSine; 
    nMC=nBoot; iFig=106; % nMC is the # MonteCarlo reps
    
    uStarTh=uThresh+zeros(length(NEE),1); % tjk, Aug 2012
    
    % nMC usually set to 1000 but will take a long time to run, tjk.
    [ruEnsemble,ruNEE,bNeg,bPos]=abNacpFcrnRandomUncertainty20100423_2 ...
        (t,NEE,uStar,PPFD,Ta,Ts,fNight,uStarTh,DirOut,cSiteYr,nMC,iFig,k,cFlux);
    bNeg(1)=bNeg(1)*k;
    bPos(1)=bPos(1)*k; % conver the units of the intercept
    
%     eval(['save ' DirOut 'randomUncertainty/randomUncertaintyFcrn_' cSiteYr ...
%         ' cSite cYr ruEnsemble ruNEE;']);
    
    if iFig>0;
        eval(['print -djpeg100 ' DirOut 'randomUncertainty/figures/randomUncertaintyFcrn_NEE' cSiteYr ';']);
    end;
    
    %	========================================================================
   
else
    ndays=length(t)/nRecsPerDay;

    uuEnsemble.NEE.Annual=NaN;
    ruEnsemble.NEE.Annual=NaN;
    uuEnsemble.NEE.Daily=nan(1,ndays);
    ruEnsemble.NEE.Daily=nan(1,ndays);
    uuEnsemble.NEE.Hourly=nan(1,ndays*nRecsPerDay);
    ruEnsemble.NEE.Hourly=nan(1,ndays*nRecsPerDay);

    uuEnsemble.GPP.Annual=NaN;
    ruEnsemble.GPP.Annual=NaN;
    uuEnsemble.GPP.Daily=nan(ndays,1);
    ruEnsemble.GPP.Daily=nan(ndays,1);
    uuEnsemble.GPP.Hourly=nan(ndays*nRecsPerDay,1);
    ruEnsemble.GPP.Hourly=nan(ndays*nRecsPerDay,1);

    uuEnsemble.RE.Annual=NaN;
    ruEnsemble.RE.Annual=NaN;
    uuEnsemble.RE.Daily=nan(ndays,1);
    ruEnsemble.RE.Daily=nan(ndays,1);
    uuEnsemble.RE.Hourly=nan(ndays*nRecsPerDay,1);
    ruEnsemble.RE.Hourly=nan(ndays*nRecsPerDay,1);

    NeeQC=NaN;
    NEEgfgC=NaN;
    bNeg=[NaN NaN];
    bPos=[NaN NaN];
end

disp(' '); disp(' ');

diary off;


%% NOTES, tjk

% - NEE_ref_joinUnc_y and NEE_ref_joinUnc_c = join uncertainty estimation of NEE_ref (_y and _c) (random + ustar filtering) calculated for each year as [sqrt(NEE_ref_randUnc^2 + ((NEE_84-NEE_16)/2)^2)]

