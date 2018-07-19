function [ruEnsemble_H, uuEnsemble_H, HQC,Hgf,bNeg,bPos] ...
    = abScriptNacpUncertaintyFcrn_fixedUstarDP_H(dataL2_4alan,cSite,cYr,columns,nBoot,uThresh)

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
cFlux='H';
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

mkdir([DirOut 'randomUncertainty/']);
% mkdir([DirOut 'uStarThAnalysis/']);
% mkdir([DirOut 'uStarTh2NeeUncertainty/']);

mkdir([DirOut 'randomUncertainty/figures/']);
% mkdir([DirOut 'uStarThAnalysis/figures/']);
% mkdir([DirOut 'uStarTh2NeeUncertainty/figures/']);

%	========================================================================

% 	cSite='UNkwn'; cYr='9999';
cSiteYr=[cSite '-' cYr];

FileDiary=[DirLog cSiteYr '.log']; diary(FileDiary);



% dataL2_4alan(dataL2_4alan==-9999)=NaN; % set missing to NaN




% 	assign
year=dataL2_4alan(:,columns.year);
ddoy=dataL2_4alan(:,columns.ddoy);

% 	Fc=dataL2_4alan(:,columns.Fc);
H=dataL2_4alan(:,columns.H);
% rough LE filter for bad data
H(H<-999)=NaN;
H(H>999)=NaN;

uStar=dataL2_4alan(:,columns.ustar);
Ta=dataL2_4alan(:,columns.Tair);
Ts=dataL2_4alan(:,columns.SoilT);
PPFD=dataL2_4alan(:,columns.PAR);
Vpd=dataL2_4alan(:,columns.VPD);

H(H==-9999)=NaN;
uStar(uStar==-9999)=NaN;
Ta(Ta==-9999)=NaN;
Ts(Ts==-9999)=NaN;
PPFD(PPFD==-9999)=NaN;
Vpd(Vpd==-9999)=NaN;


%	derive and accomodate

t=datenum(year,ones(length(year),1),ddoy); nt=length(t);
nRecsPerDay=round(1/nanmedian(diff(t)));
% if length(ddoy)>1700
%     nRecsPerDay=48;
% else
%     nRecsPerDay=24;
% end
nMinsPerDay=24*60/nRecsPerDay;

nMinsPerDay2=60*60*(24)/nRecsPerDay; 
% convert from W m-2 to mm
k=1; %nMinsPerDay2/2454000;

%(60*30)*tmpLE/2454000;
% W m-2 = J s-1
% kg m-2 s-1 = J s-1/2454000 (= j s-1/latent heat of vaporization [2.5*10^6 J kg-1])
% kg m-2 hh-1 = kg m-2 s-1*(60*30)
        

% change here to change the night time detection limit
% Ts should be Tsoil 0 need to use gap filled Tsoil, tjk
fNight=PPFD<20; % had to relax from 10 to 20
% 	Ts=Ta; % Ts has gaps; prefer to use shallow Ts or mean of Ta and Ts.

%	========================================================================

%	Mildly filter H to exclude extreme outliers

% 	nsP=7; nsCi=7; ndWindow=28; iFig=101;
nsP=5; nsCi=7; ndWindow=20; iFig=0;

[HQC,iOut] = abNacpNeeOutlierExclusion ...
    (t,H,PPFD,Ta,Ts,fNight,nsP,nsCi,ndWindow,cSiteYr,iFig,cFlux);

nOut=length(iOut);

disp(' ');
disp(sprintf([cSiteYr ' %4.0f outliers excluded.'], nOut));
disp(' ');

if iFig>0
    eval(['print -djpeg100 ' DirQC cFlux 'Outliers_' cSiteYr ';']);
end

H=HQC;
%HQC=k*HQC;

%	========================================================================

%	Plot input data for visual inspection.

iFig=102;

abNacpPlotUncertaintyInputs_2(t,H,uStar,PPFD,Ta,Ts,Vpd,fNight,cSiteYr,iFig,cFlux);

if iFig>0
    eval(['print -djpeg100 ' DirQC 'PlotInputs_H' cSiteYr ';']);
end

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

if nansum(H)~=0 && ~isnan(nansum(H)) % to deal with Howland not having any 2009 data
    FracSEBClosure=1; muStarTh=mean(uStarThs,2); iFig=110;
        
    % gap fill the data once and plot the gap-filling figure at mean ustar
    % threshold, tjk
  [NEP,NEPgf,NEPhat] = ...
        gapFill_LE_H_ANN_tkSept2013 ...
        (t,H,uStar,PPFD,Ta,Ts,PPFD,Ta,Ts,fNight,muStarTh,FracSEBClosure,cSiteYr,iFig);
      
    
    aNEP=k*nansum(NEPgf); nNEP=sum(~isnan(NEPgf)); nMissNEP=sum(isnan(NEPgf));
%     aRE=k*nansum(REgf); nMissRE=sum(isnan(REgf));
%     aGPP=k*nansum(GPPgf); nMissGPP=sum(isnan(GPPgf));
%     
%     LEgfgC=k*NEPgf;
    Hgf=k*NEPgf;
    Hhat=NEPhat;
    
    disp(sprintf('Filling gaps in %s, u*Th %5.3f   %4.0f %4.0f %4.0f   n %5.0f nMiss %4.0f %4.0f %4.0f', ...
        cSiteYr, mean(muStarTh), aNEP,nNEP, nMissNEP)); disp(' ');
    
    if iFig>0;
        eval(['print -djpeg100 ' DirQC 'GapFillingFcrn_H' cSiteYr ';']);
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
    
    
    uuEnsemble_H= [];
    uuH = [];
  
%     [uuEnsemble_H,uuH] = ...
%         abNacpFcrnUStarTh2NepUncertainty20100423 ...
%         (t,H,uStar,PPFD,Ta,Ts,PPFD,Ta,Ts,fNight,uStarThs,cSiteYr,iFig,k);
%     
%     eval(['save ' DirOut 'uStarTh2NeeUncertainty/uStarTh2NeeFcrnUncertainty_' cSiteYr ...
%         ' cSite cYr uuEnsemble_H,uuH;']);
    
    if iFig>0;
        eval(['print -djpeg100 ' DirOut 'uStarTh2NeeUncertainty/figures/uStarTh2HFcrnUncertainty_' cSiteYr ';']);
    end
    
    
    
    %	========================================================================
    
    %	Random uncertainty analysis
    
    %	Package as function separately for NEE and GPP RE
    %	have entire array output as an option.
    
    % % % 	uStarTh=mean(CpBoot);
    
    
    %uStarTh=CpSine; 
    nMC=nBoot; iFig=106; % nMC is the # MonteCarlo reps
    
    uStarTh=uThresh+zeros(length(H),1); % tjk, Aug 2012
    
    % nMC usually set to 1000 but will take a long time to run, tjk.
    
    [ruEnsemble_H,ruH,bNeg,bPos]=abNacpFcrnRandomUncertainty20100423_LE2 ...
         (t,H,Hhat,uStar,PPFD,Ta,Ts,fNight,uStarTh,DirOut,cSiteYr,nMC,iFig,k,cFlux);
%         (t,H,uStar,PPFD,Ta,Ts,fNight,uStarTh,DirOut,cSiteYr,nMC,iFig,k,cFlux);
   bNeg(1)=bNeg(1)*k;
    bPos(1)=bPos(1)*k; % conver the units of the intercept
  
%     eval(['save ' DirOut 'randomUncertainty/randomUncertaintyFcrn_' cSiteYr ...
%         ' cSite cYr ruEnsemble_H ruH;']);
     
    if iFig>0;
        eval(['print -djpeg100 ' DirOut 'randomUncertainty/figures/randomUncertaintyFcrn_H' cSiteYr ';']);
    end;
    
    %	========================================================================
   
else
    ndays=length(t)/nRecsPerDay;

    uuEnsemble_H.NEE.Annual=NaN;
    ruEnsemble_H.NEE.Annual=NaN;
    uuEnsemble_H.NEE.Daily=nan(ndays,1);
    ruEnsemble_H.NEE.Daily=nan(ndays,1);
    uuEnsemble_H.NEE.Hourly=nan(ndays*nRecsPerDay,1);
    ruEnsemble_H.NEE.Hourly=nan(ndays*nRecsPerDay,1);
    
    HQC=NaN;
    Hgf=NaN;
    bNeg=[NaN NaN];
    bPos=[NaN NaN];
end

disp(' '); disp(' ');

diary off;


%% NOTES, tjk

% - NEE_ref_joinUnc_y and NEE_ref_joinUnc_c = join uncertainty estimation of NEE_ref (_y and _c) (random + ustar filtering) calculated for each year as [sqrt(NEE_ref_randUnc^2 + ((NEE_84-NEE_16)/2)^2)]

