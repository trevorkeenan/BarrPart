function [NEP,R,GEP,NEPgf,Rgf,GEPgf,RHat,GEPHat,RHat0,GEPHat0,mtR,mcR,mtGEP,mcGEP,bR,bGEP] = ...
    abNacpFcrnCO2Flux2NEP20090205 ...
    (t,NEE,uStar,PPFD,Ta,Ts,PPFDGF,TaGF,TsGF,isNight,uStarTh,FracSEBClosure,cSiteYr,iFig);

%abNacpFcrnCO2Flux2NEP20090205
%
%created 29 Sept 2010 to reset to normal FCRN program.
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
% 1. NEP is estimated as -NEE after excluding low u* data at night and
% 	applying the optional energy closure adjustment.
% 2. Measured R is estimated from NEE at night and during the non-GrowingSeason
% 3. An empirical R=f(Ts,t) model is fit based on the data from 2
% 4. The R=f(Ts,t) model from 3 is used to estimate R during the day
% 	and to fill gaps in R at night
% 5. GEP is estimated as NEP+R (daytime, growing season) or zero
% 	(nighttime and non-GrowingSeason)
% 6. An empirical GEP=f(PPFD,t) model is fit to the GEP data from 5
% 7. The model from 6 is used to fill gaps in daytime, growing-season GEP
% 8. Gaps in NEP are filled using modelled GEP-R
%
% 	A moving window is used to capture the time variation in
% 	the R and GEP model parameters cR(t) and cGEP(t):
% 		R 		= cR(t)*bR(1)/(1+exp[bR(2)*(bR(3)-Ts)]);
% 		GEP 	= cGEP(t)*bGEP(1)*PPFD/(bGEP(2)+PPFD);
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

%	Several non-standard MatLab functions are celled:
%	- doy
%	- fcFillSmallGapsByLinInterp
%	- fr_function_fit
%	- gapsize
%	- fcdatevec
%	- fcnanmean
%	- fcnanmedian
%	- fcnansum
%	- fcNanInterp1

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

%	Fill small gaps in input data.
%	NEP done AFTER low-u* rejection.

warning off;

%GapSizeMax = 4;
GapSizeMax = 5;

uStar = fcFillSmallGapsByLinInterp(uStar,GapSizeMax);
PPFD = fcFillSmallGapsByLinInterp(PPFD,GapSizeMax);
Ta = fcFillSmallGapsByLinInterp(Ta,GapSizeMax);
Ts = fcFillSmallGapsByLinInterp(Ts,GapSizeMax);
PPFDGF = fcFillSmallGapsByLinInterp(PPFDGF,GapSizeMax);
TaGF = fcFillSmallGapsByLinInterp(TaGF,GapSizeMax);
TsGF = fcFillSmallGapsByLinInterp(TsGF,GapSizeMax);

%	========================================================================

%	Initial Assignments.

%	Assign the moving window parameters.

nRecsPerDay=round(1/nanmedian(diff(t)));

WindowWidthR = 100; WindowIncR = 20;
WindowWidthGEP = 100; WindowIncGEP = 20;

if nRecsPerDay==24;
    WindowWidthR = WindowWidthR/2; WindowIncR = WindowIncR/2;
    WindowWidthGEP = WindowWidthGEP/2; WindowIncGEP = WindowIncGEP/2;
end;

%	Assign conditions (is*) and indices (i*).

nt = length(NEE);

isDay = ~ isNight;
coldA=mean(TaGF)-2*std(TaGF);
coldS=mean(TsGF)-2*std(TsGF);
isColdSeason = TaGF<=coldA & TsGF<=coldS;
isWarmSeason = TaGF>coldA & TsGF>coldS;
isNotColdSeason = ~ isColdSeason;

iNight = find(isNight);
iDay = find(isDay);
iColdSeason = find(isColdSeason);
iWarmSeason = find(isWarmSeason);
iNotColdSeason = find(isNotColdSeason);
iLowuStar = find(uStar<uStarTh);

%	========================================================================

%	Assign NEP.

NEP = -NEE;

%	Exclude low-u* data at night.

iReject = intersect(iLowuStar,iNight);
NEP(iReject) = NaN;

%	Apply SEB-closure adjustment to the data.

NEP = NEP./FracSEBClosure;

% fixed error where filled small gaps were assigned to NEP not NEPgf.
% carried through to Rgf and GEPgf.

NEPgf=NEP; GapSizeMax = 4; NEPgf = fcFillSmallGapsByLinInterp(NEPgf,GapSizeMax);

%	========================================================================

%	Estimate and gap-fill R.

MT = NaN*ones(size(NEP));
R = MT; RHat = MT; Rgf=MT;

%	Assign R = NEP at night and during daytime periods in the cold season.

R(iNight) = -NEP(iNight);
R(iColdSeason) = -NEP(iColdSeason);
Rgf(iNight) = -NEPgf(iNight);
Rgf(iColdSeason) = -NEPgf(iColdSeason);

% R(R<-1)=NaN;  % commented out as was causing major problems by biasing
% the distribution of night time flux high. This resulted in Reco being
% overestimated by about 30%
% Rgf(Rgf<-1)=NaN;    % remove extreme outliers in R

% % % 	plot(t,Rgf,'o'); fcDateTick(t,'Mo',4,1); pause;

%	Tr used to fit R=f(Tr) is set to the mean of Ts and Ta.

Tr=(Ta+Ts)/2;
TrGF=(TaGF+TsGF)/2;

%	Fit the R model R = cR(t)*bR(1)/(1+exp[bR(2)*(bR(3)-Tr)]).

%	First, fit the bR parameters using binned means of all of the data.

iReg = find(~isnan(R+TrGF));
nReg = length(iReg);
bRGuess = [1 1 10];
lasterror('reset');
bR = nlinfit(TrGF(iReg),R(iReg),@fcCalcTs2RLogistic,bRGuess);
err = lasterror;
if ~isempty(err.message);
    disp(' '); disp('*** fcCalcTs2RLogistic failed. ***')
    disp(err.message'); disp(' ');
end; % linear regression failed.

RHat0 = fcCalcTs2RLogistic(bR,TrGF);

%	Second, use a moving window to estimate mcR(t), based on measured R
%	and the preliminary estimate of modelled R (RHat0). The values (mcR,mtR)
%	are first fit to each window and then interpolated to (cR(t),t).

%	limit data to close to the mode

% ndSwath=31; % limit to points within ndSwath days of the time mode.
% nRegN=30; % if resulting # points are less then 30 skip.

% RELAXED FOR TK
	ndSwath=5*31; % limit to points within ndSwath days of the time mode.
	nRegN=15; % if resulting # points are less then 20 skip.

i = 0; mtR = []; mcR = [];
for jReg1 = 1:WindowIncR:nReg;
    jReg2 = min(jReg1+WindowWidthR-1,nReg);
    if jReg2<=nReg;
        it=iReg(jReg1:jReg2);
        % restrict to point ndSwath days wihin the mode.
        iOut=find(abs(t(it)-mode(t(it)))>ndSwath);
        it(iOut)=[]; nit=length(it);
        if nit>=nRegN;
            i=i+1;
            mtR(i) = mean(t(it));
            mcR(i) = mean(R(it))/mean(RHat0(it));
        end;
    end; % if jReg2< = nReg;
end; % for jReg1 = 1:WindowIncR:nReg;

mcR(mcR<0)=0;


% quality control for badly fit data
med=nanmedian(mcR);
if isnan(med)
    med=median(mcR); % not sure why but nanmedian returns nan sometimes. ?$%(^??
end

stdev=25;

mcR(abs(mcR)>med+stdev)=med;

mmcR=nanmedian(mcR); 
if isnan(mmcR)
    mmcR=median(mcR); 
end    
qmcR=fcNanIqr(fcx2colvec(mcR));
% % % 	clf; hist(mcR); pause;
ns=20; iEx=find(abs(mcR-mmcR)>ns*qmcR); mtR(iEx)=[]; mcR(iEx)=[];

%	Interpolate (mtR,mcR) coefficients to (t,cR) 30 min periods.

cR = fcNanInterp1(mtR,mcR,t,'linear');

%	Set large gaps to zero.

%	If the data begin or end part-way through the time series,
%	set cR to 1.0 at the beginning and/or end. Allow a 30-day buffer.

% fix for bad fitting at Santa_Rita_Mesquite_Savanna
if strcmp(cSiteYr,'Santa_Rita_Mesquite_Savanna-2006') ||...
   strcmp(cSiteYr,'Santa_Rita_Mesquite_Savanna-2008')
    cR(1:4000)=0.5;
    cR(14000:end)=0.5;
end


%	Eliminate large gaps or gaps at the start or end of the time series.

nGapSize=fcGapSize(R);
iEx=find(nGapSize>ndSwath*nRecsPerDay); cR(iEx)=NaN;
itN = min(iReg); if itN>30*nRecsPerDay; it = 1:(itN-1); cR(it) = NaN; end;
itX = max(iReg); if itX<(nt-30*nRecsPerDay); it = (itX+1):nt; cR(it) = NaN; end;

%	Calculate modeled R (RHat).

RHat 					= cR.*RHat0;
RHat(isnan(RHat))=RHat0(isnan(RHat));
% % % 	plot(mtR,mcR,'o-'); fcDateTick(t,'Mo',4,1); pause

%	Merge R and RHat to produce Rgf (gap-filled).

igf = find(isnan(Rgf) & ~isnan(RHat));
Rgf(igf) = RHat(igf);

%	========================================================================

%	Estimate and gap-fill GEP.

%	Estimate measured GEP. Convention: NEP = GEP-R so GEP = NEP+R.

GEP = NEP+Rgf;
GEPgf = NEPgf+Rgf;

%	Set GEP to zero at night and during the cold season.

GEP(iNight) = 0; GEP(iColdSeason) = 0;
GEPgf(iNight) = 0; GEPgf(iColdSeason) = 0;

%	Fit the GEP model in two steps.

%	First, estimate the constant bGEP parameters using warm-season data only.

iReg = find(isDay & isWarmSeason & ~isnan(PPFD+GEP));
bGEPGuess = [20 100];
lasterror('reset');
tmp=GEP(iReg);
tmp(~isfinite(tmp))=NaN;
bGEP = nlinfit(PPFD(iReg),tmp,@fcCalcRs2GEP,bGEPGuess);
err = lasterror;
if ~isempty(err.message);
    disp(' '); disp('*** fcCalcRs2GEP failed. ***')
    disp(err.message'); disp(' ');
end; % linear regression failed.
GEPHat0 = fcCalcRs2GEP(bGEP,PPFDGF);

%	Second, estimate mcGEP(t) for the whole annual cycle based on the NotColdSeason GEP data.
%	Using a moving window, estimate mcGEP(t) by forced origin linear regression
% 	of measured GEP versus GEPHat0, i.e., GEP(t)  =  mcGEP(t) * GEPHat0.

iReg = find(isDay & isNotColdSeason & ~isnan(PPFD+GEP));
nReg = length(iReg);

i=0; mtGEP = []; mcGEP = [];
for jReg1 = 1:WindowIncGEP:nReg;
    jReg2 = min(jReg1+WindowWidthGEP-1,nReg);
    if jReg2<=nReg;
        it = iReg(jReg1:jReg2);
        % restrict to point ndSwath days wihin the mode.
        iOut=find(abs(t(it)-mode(t(it)))>ndSwath);
        it(iOut)=[]; nit=length(it);
        if nit>=nRegN;
            i=i+1;
            mtGEP(i) = mean(t(it));
            mcGEP(i) = mean(GEP(it))/mean(GEPHat0(it));
        end;
    end; % if jReg2< = nReg;
end; % for jReg1
mcGEP(mcGEP<0)=0;

mmcGEP=nanmedian(mcGEP); qmcGEP=fcNanIqr(fcx2colvec(mcGEP));
ns=20; iEx=find(abs(mcGEP-mmcGEP)>ns*qmcGEP); mtGEP(iEx)=[]; mcGEP(iEx)=[];

%	Interpolate (mtGEP,mcGEP) to (t,cGEP) and calculate modeled GEP (GEPHat).
%	First, drop missing values of mcGEP.

iNaN = find(isnan(mcGEP)); mtGEP(iNaN) = []; mcGEP(iNaN) = [];
cGEP = fcNanInterp1(mtGEP,mcGEP,t,'linear');

%	Eliminate large gaps or gaps at the start or end of the time series.
%	fZeroGEP causes the gapsize to be based on daytime periods only.  Added
%	5 August 2009.

fZero=ones(size(GEP)); fZero(iNight)=NaN; fZero(iColdSeason)=NaN;

nGapSize=fcGapSize(fZero.*GEP);
iEx=find(nGapSize>ndSwath*nRecsPerDay); cGEP(iEx)=NaN;
%	itN = min(iReg); if itN>30*nRecsPerDay; it = 1:(itN-1); cGEP(it) = NaN; end;
%	itX = max(iReg); if itX<(nt-30*nRecsPerDay); it = (itX+1):nt; cGEP(it) = NaN; end;

%	Also set cGEP to zero for periods when GEP is zero.
%	Probably redundant because GEP is already set to zero for these periods.

iZero = union(iNight,iColdSeason);
cGEP(iZero) = 0;

%	Estimate modelled GEP (GEPHat).

GEPHat = cGEP.*GEPHat0;
GEPHat(isnan(GEPHat))=GEPHat0(isnan(GEPHat)); % cGEP gaps issue at end of year from fcNanInterp1

%	Merge GEP and GEPHat to create GEPgf (gap-filled).

igf = find(isnan(GEPgf) & ~isnan(GEPHat));
GEPgf(igf) = GEPHat(igf);

%	========================================================================

%	Assign final gap-filled NEP = GEP-R values and fill gaps with modeled estimates.

igf = find(isnan(NEPgf) & ~isnan(Rgf+GEPgf));
NEPgf(igf) = GEPgf(igf)-Rgf(igf);

% check for massive outliers
NEPgf(abs(NEPgf)>9999)=0;
Rgf(abs(Rgf)>9999)=0;
GEPgf(abs(GEPgf)>9999)=0;


%	========================================================================
%	========================================================================

if iFig>0;
    
    fcFigLoc(iFig,0.6,0.9,'SE');
    
    subplot('position',[0.05 0.76 0.42 0.19]); hold on;
    plot(t,-NEPgf,'r.',t,-NEP,'b.'); fcDateTick(t,'Mo',4,1); xlim([min(t) max(t)]);
    title(sprintf('%s  NEE nMiss %g',cSiteYr,sum(isnan(fcx2rowvec(NEPgf)))));
    subplot('position',[0.05 0.52 0.42 0.19]);
    plot(t,Rgf,'r.',t,R,'b.'); fcDateTick(t,'Mo',4,1); xlim([min(t) max(t)]);
    title(sprintf('RE nMiss %g',sum(isnan(fcx2rowvec(Rgf)))));
    subplot('position',[0.05 0.28 0.42 0.19]);
    plot(t,GEPgf,'r.',t,GEP,'b.'); fcDateTick(t,'Mo',4,1); xlim([min(t) max(t)]);
    title(sprintf('GPP nMiss %g',sum(isnan(fcx2rowvec(GEPgf)))));
    subplot('position',[0.05 0.04 0.42 0.19]);
    plot(t,[TaGF TsGF],'.'); fcDateTick(t,'Mo',4,1); xlim([min(t) max(t)]);
    title(sprintf('Ta Ts nMiss %g',sum(isnan([TaGF TsGF]))));
    
    subplot('position',[0.55 0.76 0.42 0.19]); hold on;
    plot(TrGF,Rgf,'r.',TrGF,R,'b.');
    hold on; plot(min(TrGF):0.01:max(TrGF),fcCalcTs2RLogistic(bR,min(TrGF):0.01:max(TrGF)),'c-');
    xlim([min(TrGF) max(TrGF)]);
    ylim([nanmin(R) nanmax(R)]);
    box on; grid on;
    title(sprintf('RE vs T,  b = %3.1f  %3.2f  %3.1f',bR));
    
    subplot('position',[0.55 0.52 0.42 0.19]);
    tmpPPFDGF=PPFDGF(iWarmSeason);
    tmpGEP=GEP(iWarmSeason);
    tmpGEPgf=GEPgf(iWarmSeason);
    
%    figure;plot(PPFDGF(iDay),GEPgf(iDay),'r.',PPFDGF(iDay),GEP(iDay),'b.');
    plot(tmpPPFDGF,tmpGEPgf,'r.',tmpPPFDGF,tmpGEP,'b.');
    hold on; plot(0:max(PPFDGF),fcCalcRs2GEP(bGEP,0:max(PPFDGF)),'c-');
    xlim([min(PPFDGF) max(PPFDGF)]);
    ylim([min(GEP) max(GEP)]);
    box on; grid on;
    title(sprintf('GPP vs PPFD,  b = %3.1f  %3.1f  %3.1f',bGEP));
    
    subplot('position',[0.55 0.28 0.42 0.19]);
    plot(mtR,mcR,'o'); fcDateTick(mtR,'Mo',4,1); xlim([min(mtR) max(mtR)]);
    title('cR'); ylim([-1 3]);
    subplot('position',[0.55 0.04 0.42 0.19]);
    plot(mtGEP,mcGEP,'o'); fcDateTick(mtGEP,'Mo',4,1); xlim([min(mtGEP) max(mtGEP)]);
    title('cGPP');
    
end;

%	========================================================================
%	========================================================================

