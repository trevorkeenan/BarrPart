function [ruEnsemble_LE, uuEnsemble_LE, LEQC,LEgfgC,bNeg,bPos] ...
    = abScriptNacpUncertaintyFcrn_fixedUstarDP_LE(dataL2_4alan,cSite,cYr,columns,nBoot,uThresh)

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
cFlux='LE';

close all;
%     clear all;

warning off;

% addpath(genpath('./functions/'),'./functions2/');


mkdir('../../Folder-Structure_Results/Figures/'); 
DirRoot='../../Folder-Structure_Results/Figures/';
DirRoot2='../../Folder-Structure_Results/Figures/';

mkdir([DirRoot2 'Output/FCRN/']); DirOut=[DirRoot2 'Output/FCRN/'];
mkdir([DirRoot2 'figuresQC/']); DirQC=[DirRoot2 'figuresQC/'];
mkdir([DirRoot2 'log/']); DirLog=[DirRoot2 'Log/'];

mkdir([DirOut 'randomUncertainty/']);

mkdir([DirOut 'randomUncertainty/figures/']);

%	========================================================================

cSiteYr=[cSite '-' cYr];

FileDiary=[DirLog cSiteYr '.log']; diary(FileDiary);


% 	assign
year=dataL2_4alan(:,columns.year);
ddoy=dataL2_4alan(:,columns.ddoy);

LE=dataL2_4alan(:,columns.LE);

% rough LE filter for bad data
LE(LE<-100)=NaN;
LE(LE>9999)=NaN;

uStar=dataL2_4alan(:,columns.ustar);
Ta=dataL2_4alan(:,columns.Tair);
Ts=dataL2_4alan(:,columns.SoilT);
PPFD=dataL2_4alan(:,columns.PAR);
Vpd=dataL2_4alan(:,columns.VPD);

LE(LE==-9999)=NaN;
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
k=nMinsPerDay2/2454000;

%(60*30)*tmpLE/2454000;
% W m-2 = J s-1
% kg m-2 s-1 = J s-1/2454000 (= j s-1/latent heat of vaporization [2.5*10^6 J kg-1])
% kg m-2 hh-1 = kg m-2 s-1*(60*30)
        

% change here to change the night time detection limit
% Ts should be Tsoil 0 need to use gap filled Tsoil, tjk
fNight=PPFD<20; % had to relax from 10 to 20
% 	Ts=Ta; % Ts has gaps; prefer to use shallow Ts or mean of Ta and Ts.

%	========================================================================

%	Mildly filter LE to exclude extreme outliers
nsP=5; nsCi=7; ndWindow=20; iFig=0;

[LEQC,iOut] = abNacpNeeOutlierExclusion ...
    (t,LE,PPFD,Ta,Ts,fNight,nsP,nsCi,ndWindow,cSiteYr,iFig,cFlux);
nOut=length(iOut);

disp(' ');
disp(sprintf([cSiteYr ' %4.0f outliers excluded.'], nOut));
disp(' ');

if iFig>0;
    eval(['print -djpeg100 ' DirQC cFlux 'Outliers_' cSiteYr ';']);
end;

% LEQC=LE;

LE=LEQC;
LEQC=k*LEQC;

%	========================================================================

%	Plot input data for visual inspection.

iFig=102;

abNacpPlotUncertaintyInputs_2(t,LE,uStar,PPFD,Ta,Ts,Vpd,fNight,cSiteYr,iFig,cFlux);

if iFig>0;
    eval(['print -djpeg100 ' DirQC 'PlotInputs_Le' cSiteYr ';']);
end;

%	========================================================================


% tjk, Aug 2012
CpBoot=uThresh;
Cp=uThresh;
    
    %	========================================================================
    
    %	check processing at mean uStarTh.
    
    
uStarThs = Cp;  % tjk, Aug 2012

if nansum(LE)~=0 && ~isnan(nansum(LE)) % to deal with Howland not having any 2009 data
    FracSEBClosure=1; muStarTh=mean(uStarThs,2); iFig=109;
        
    % gap fill the data once and plot the gap-filling figure at mean ustar
    % threshold, tjk
%     [NEP,RE,GPP,NEPgf,REgf,GPPgf,RHat,GPPHat,RHat0,GPPHat0,mtR,mcR,mtGPP,mcGPP,bR,bGPP] = ...
%         abNacpFcrnCO2Flux2NEP20090205 ...
%         (t,LE,uStar,PPFD,Ta,Ts,PPFD,Ta,Ts,fNight,muStarTh,FracSEBClosure,cSiteYr,iFig);
    [NEP,NEPgf,NEPhat] = ...
        gapFill_LE_H_ANN_tkSept2013 ...
        (t,LE,uStar,PPFD,Ta,Ts,PPFD,Ta,Ts,fNight,muStarTh,FracSEBClosure,cSiteYr,iFig);
    
    
    aNEP=k*nansum(NEPgf); nNEP=sum(~isnan(NEPgf)); nMissNEP=sum(isnan(NEPgf));
%     aRE=k*nansum(REgf); nMissRE=sum(isnan(REgf));
%     aGPP=k*nansum(GPPgf); nMissGPP=sum(isnan(GPPgf));
    
    LEgfgC=k*NEPgf;
    LEhat=NEPhat;
    
    disp(sprintf('Filling gaps in %s, u*Th %5.3f   %4.0f %4.0f %4.0f   n %5.0f nMiss %4.0f %4.0f %4.0f', ...
        cSiteYr, mean(muStarTh), aNEP, nNEP, nMissNEP)); disp(' ');
    
    if iFig>0;
        eval(['print -djpeg100 ' DirQC 'GapFillingFcrn_LE' cSiteYr ';']);
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
    
    
    
    uuEnsemble_LE = [];
    uuLE = []; %...
%         abNacpFcrnUStarTh2NepUncertainty20100423 ...
%         (t,LE,uStar,PPFD,Ta,Ts,PPFD,Ta,Ts,fNight,uStarThs,cSiteYr,iFig,k);
        
    if iFig>0;
        eval(['print -djpeg100 ' DirOut 'uStarTh2NeeUncertainty/figures/uStarTh2LEFcrnUncertainty_' cSiteYr ';']);
    end
    
    
    
    %	========================================================================
    
    %	Random uncertainty analysis
    
    %	Package as function separately for NEE and GPP RE
    %	have entire array output as an option.
    
    nMC=nBoot; iFig=106; % nMC is the # MonteCarlo reps
    
    uStarTh=uThresh+zeros(length(LE),1); % tjk, Aug 2012
    
    % nMC usually set to 1000 but will take a long time to run, tjk.
    
    [ruEnsemble_LE,ruLE,bNeg,bPos]=abNacpFcrnRandomUncertainty20100423_LE2 ...
        (t,LE,LEhat,uStar,PPFD,Ta,Ts,fNight,uStarTh,DirOut,cSiteYr,nMC,iFig,k,cFlux);
    bNeg(1)=bNeg(1)*k;
    bPos(1)=bPos(1)*k; % conver the units of the intercept
  
     
    if iFig>0;
        eval(['print -djpeg100 ' DirOut 'randomUncertainty/figures/randomUncertaintyFcrn_LE' cSiteYr ';']);
    end;
    
    %	========================================================================
   
else
    ndays=length(t)/nRecsPerDay;

    uuEnsemble_LE.NEE.Annual=NaN;
    ruEnsemble_LE.NEE.Annual=NaN;
    uuEnsemble_LE.NEE.Daily=nan(ndays,1);
    ruEnsemble_LE.NEE.Daily=nan(ndays,1);
    uuEnsemble_LE.NEE.Hourly=nan(ndays*nRecsPerDay,1);
    ruEnsemble_LE.NEE.Hourly=nan(ndays*nRecsPerDay,1);

    LEQC=NaN;
    LEgfgC=NaN;
    bNeg=[NaN NaN];
    bPos=[NaN NaN];
end

disp(' '); disp(' ');

diary off;
