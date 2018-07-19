function [FLUX,FLUXgf,FLUXhat] = ...
    gapFill_LE_H_ANN_tkSept2013 ...
    (t,FLUX,uStar,PPFD,Ta,Ts,PPFDGF,TaGF,TsGF,isNight,uStarTh,FracSEBClosure,cSiteYr,iFig,cFlux);

% script based on original code in: abNacpFcrnCO2Flux2FLUX20090205
% modified to use ANN instead of RE and GPP partitioning
% for use with LE and H

% 1. Gap fill night time LE or H fluxes using soil temperature
% 2. Gap fill day time LE and H fluxes using air temperature, PPFD, soil T.


%Syntax:
%------
%
%[FLUX,FLUXgf]
%	= CO2Flux2FLUX(t,FLUX,uStar,PPFD,Ta,Ts,PPFDGF,TaGF,TsGF,isNight,uStarTh,FracSEBClosure);
%
%Input Arguments:
%---------------
%
%-	t is the decimal day vector
%-	flux is the measured flux (LE or H) vector
%-	uStar is the friction velocity vector
%-	PPFD is the downwelling PAR flux density vector (used to model GEP)
%-	Ta is the air temperature vector (used optionally to model GEP)
%-	Ts is the soil temperature vector (used to model R)
%-	PPFDGF, TaGF, and TsGF are gap-filled values of PPFD, Ta, and Ts
%-	isNight is a vector that indicates whether itReg is day (0) or night (1)
%-	uStarTh is the scalar (or vector) uStar threshold below which nighttime FLUX data are rejected
%-	FracSEBClosure is a scalar or vector SEB-closure adjustment divider
%-	LogFile is the (optional) log file directory and name; create if not empty.
%-	Plots (optional) specifies the plots to create.

%	========================================================================
%	========================================================================

%	Notes

%	Written by Trevor Keenan (September 2013)


%	========================================================================
%	========================================================================


%	Fill small gaps in input data.
%	FLUX done AFTER low-u* rejection.

warning off;

GapSizeMax = 4;

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

% define an X day window
WindowInc = 180*nRecsPerDay;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Define a Neural Network
% Create a Fitting Network
            hiddenLayerSize = 5;
            
            net1 = fitnet(hiddenLayerSize);
            
            % Choose Input and Output Pre/Post-Processing Functions
            % For a list of all processing functions type: help nnprocess
            net1.inputs{1}.processFcns = {'removeconstantrows','mapminmax'};
            net1.outputs{2}.processFcns = {'removeconstantrows','mapminmax'};
            
            % Setup Division of Data for Training, Validation, Testing
            % For a list of all data division functions type: help nndivide
            net1.divideFcn = 'dividerand';  % Divide data randomly
            net1.divideMode = 'sample';  % Divide up every sample
            
            training=60; % 70
            validating=20; % 15
            testing=20; % 15
            
            net1.divideParam.trainRatio = training/100;
            net1.divideParam.valRatio = validating/100;
            net1.divideParam.testRatio = testing/100;
            
            % For help on training function 'trainlm' type: help trainlm
            % For a list of all training functions type: help nntrain
            net1.trainFcn = 'trainlm';  % Levenberg-Marquardt
            
            % Choose a Performance Function
            % For a list of all performance functions type: help nnperformance
            net1.performFcn = 'mse';  % Mean squared error
            
            % Choose Plot Functions
            % For a list of all plot functions type: help nnplot
            net1.plotFcns = {'plotperform','plottrainstate','ploterrhist', ...
                'plotregression', 'plotfit'};
            
            net1.trainParam.epochs = 50;
            net1.trainParam.goal = 0.0001;
            
            net1.trainParam.showWindow = false;
            net1.trainParam.showCommandLine = false;
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%	Assign conditions (is*) and indices (i*).

nt = length(FLUX);

isDay = ~ isNight;
isColdSeason = TaGF<=0 & TsGF<=0;
isWarmSeason = TaGF>0 & TsGF>0;
isNotColdSeason = ~ isColdSeason;

iNight = find(isNight);
iDay = find(isDay);
iColdSeason = find(isColdSeason);
iWarmSeason = find(isWarmSeason);
iNotColdSeason = find(isNotColdSeason);
iLowuStar = find(uStar<uStarTh);

%	========================================================================

%	Exclude low-u* data at night.

iReject = intersect(iLowuStar,iNight);
FLUX(iReject) = NaN;

FLUX(FLUX<-175)=NaN;    % hard filter for very large negative outliers

%	Apply SEB-closure adjustment to the data.
FLUX = FLUX./FracSEBClosure;

% fixed error where filled small gaps were assigned to FLUX not FLUXgf.

FLUXgf=FLUX; GapSizeMax = 4; FLUXgf = fcFillSmallGapsByLinInterp(FLUXgf,GapSizeMax);

% do not use negative FLUX values for gap filling
FLUX2=FLUX;
FLUX2(FLUX2<-50)=NaN;
%	========================================================================


%	Use a moving window to estimate mcR(t), based on measured R
%	and the preliminary estimate of modelled R (RHat0). The values (mcR,mtR)
%	are first fit to each window and then interpolated to (cR(t),t).

%	limit data to close to the mode

% RELAXED FOR TK
ndSwath=5*31; % limit to points within ndSwath days of the time mode.
nRegN=15; % if resulting # points are less then 20 skip.

iReg=1:length(FLUX);
nReg = length(iReg);

% first fit the whole year night-time data
FLUXhat=FLUX+NaN;
% FLUXhat(iNight) = FLUX(iNight);

% fit night time data all year together using soil temperature
inputs = [iReg(iNight)'/length(FLUX) TsGF(iNight) TaGF(iNight)]';
% convert time to cosine
inputs(1,:)=cos(2*pi*inputs(1,:));

% define the target
targets = FLUX2(iNight)';




% Train the Network
            clear outputs1 testPerformance
            numIt=3;
            testPerformance=zeros(1,numIt);
            outputs1=zeros(numIt,length(targets));
            for tjk=1:numIt
                [net,tr] = train(net1,inputs,targets);
                outputs = net(inputs);
                
                testTargets = targets  .* tr.testMask{1};
                testPerformance(tjk) = perform(net,testTargets,outputs);
                
                outputs1(tjk,:) = outputs;
                clear outputs testTargets
            end
            
            % use output from the best performing network
            
            best=find(testPerformance==min(testPerformance));
            if length(best)>1
                % to account for the possiblility that there is more than one,
                % or no, best ANN
                best=best(1);
            end
            if length(outputs1(best,:))==length(FLUXhat(iNight))
                FLUXhat(iNight)=outputs1(best,:);
            end



% Now fit to the day-time data using air temperature
% fit night time data all year together using soil temperature
inputs = [iReg(iDay)'/length(FLUX) TaGF(iDay) TsGF(iDay) PPFDGF(iDay)]';
% convert time to cosine
inputs(1,:)=cos(2*pi*inputs(1,:));

% define the target
targets = FLUX2(iDay)';

FLUXhat(FLUXhat<-50)=0;


% Train the Network
            clear outputs1 testPerformance
            numIt=3;
            testPerformance=zeros(1,numIt);
            outputs1=zeros(numIt,length(targets));
            for tjk=1:numIt
                [net,tr] = train(net1,inputs,targets);
                outputs = net(inputs);
                
                testTargets = targets  .* tr.testMask{1};
                testPerformance(tjk) = perform(net,testTargets,outputs);
                
                outputs1(tjk,:) = outputs;
                clear outputs testTargets
            end
            
            % use output from the best performing network
            
            best=find(testPerformance==min(testPerformance));
            if length(best)>1
                % to account for the possiblility that there is more than one,
                % or no, best ANN
                best=best(1);
            end
            if length(outputs1(best,:))==length(FLUXhat(iDay))
                FLUXhat(iDay)=outputs1(best,:);
            end




% %	Fit the ANN model to different time slices in the year.
% % FLUXhat2=FLUX+NaN;
% for jReg1 = 1:WindowInc:nReg;
%     jReg2 = min(jReg1+WindowInc-1,nReg);
%     if jReg2<=nReg; 
%         it = iReg(jReg1:jReg2);
%         % restrict to point ndSwath days wihin the mode.
%         iOut=find(abs(t(it)-mode(t(it)))>ndSwath);
%         it(iOut)=[]; 
%         it(ismember(it,intersect(it,iNight)))=[];
%         nit=length(it);
%         if nit>=nRegN;
%             % put the ANN here.
%             % use AirT, PAR for training
%             inputs = [iReg(it)'/length(FLUX) PPFDGF(it) TaGF(it)]';
%             % convert time to cosine
%             inputs(1,:)=cos(2*pi*inputs(1,:));
%             
%             % define the target
%             targets = FLUX2(it)';
%             
%             
%             % Train the Network
%             clear outputs1 testPerformance
%             numIt=5;
%             testPerformance=zeros(1,numIt);
%             outputs1=zeros(numIt,length(targets));
%             for tjk=1:numIt
%                 [net,tr] = train(net1,inputs,targets);
%                 outputs = net(inputs);
%                 
%                 %                 errors = gsubtract(targets,outputs);
%                 %                 performance = perform(net,targets,outputs);
%                 %
%                 % Recalculate Training, Validation and Test Performance
%                 %                 trainTargets = targets .* tr.trainMask{1};
%                 %                 valTargets = targets  .* tr.valMask{1};
%                 testTargets = targets  .* tr.testMask{1};
%                 %                 trainPerformance = perform(net,trainTargets,outputs);
%                 %                 valPerformance = perform(net,valTargets,outputs);
%                 testPerformance(tjk) = perform(net,testTargets,outputs);
%                 
%                 outputs1(tjk,:) = outputs;
%                 clear outputs testTargets
%             end
%             
%             % use output from the best performing network
%             
%             best=find(testPerformance==min(testPerformance));
%             if length(best)>1
%                 % to account for the possiblility that there is more than one,
%                 % or no, best ANN
%                 best=best(1);
%             end
%             if length(outputs1(best,:))==length(FLUXhat(it))
%                 FLUXhat(it)=outputs1(best,:);
%             end
%             
%         end;
%     end; % if jReg2< = nReg;
% end; % for jReg1



%	========================================================================

%	Assign final gap-filled FLUX = GEP-R values and fill gaps with modeled estimates.

igf = find(isnan(FLUX));
FLUXgf(igf) = FLUXhat(igf);

%	========================================================================
%	========================================================================

if iFig>0;
%%    
    fcFigLoc(iFig,0.6,0.9,'SE');
    
    %     subplot('position',[0.05 0.76 0.42 0.19]); hold on;
    plot(t,FLUXgf,'r.',t,FLUX,'b.'); fcDateTick(t,'Mo',4,1); xlim([min(t) max(t)]);
    title(sprintf('%s  FLUX nMiss %g',cSiteYr,sum(isnan(fcx2rowvec(FLUXgf)))));
    %     subplot('position',[0.05 0.52 0.42 0.19]);
    %     plot(t,Rgf,'r.',t,R,'b.'); fcDateTick(t,'Mo',4,1); xlim([min(t) max(t)]);
    %     title(sprintf('RE nMiss %g',sum(isnan(fcx2rowvec(Rgf)))));
    %     subplot('position',[0.05 0.28 0.42 0.19]);
    %     plot(t,GEPgf,'r.',t,GEP,'b.'); fcDateTick(t,'Mo',4,1); xlim([min(t) max(t)]);
    %     title(sprintf('GPP nMiss %g',sum(isnan(fcx2rowvec(GEPgf)))));
    %     subplot('position',[0.05 0.04 0.42 0.19]);
    %     plot(t,[TaGF TsGF],'.'); fcDateTick(t,'Mo',4,1); xlim([min(t) max(t)]);
    %     title(sprintf('Ta Ts nMiss %g',sum(isnan([TaGF TsGF]))));
    
    %     subplot('position',[0.55 0.76 0.42 0.19]); hold on;
    %     plot(TrGF,Rgf,'r.',TrGF,R,'b.');
    %     hold on; plot(min(TrGF):0.01:max(TrGF),fcCalcTs2RLogistic(bR,min(TrGF):0.01:max(TrGF)),'c-');
    %     xlim([min(TrGF) max(TrGF)]);
    %     ylim([min(R) max(R)]);
    box on; grid on;
    %     title(sprintf('RE vs T,  b = %3.1f  %3.2f  %3.1f',bR));
    %     subplot('position',[0.55 0.52 0.42 0.19]);
    %     plot(PPFDGF,GEPgf,'r.',PPFDGF,GEP,'b.');
    %     hold on; plot(0:max(PPFDGF),fcCalcRs2GEP(bGEP,0:max(PPFDGF)),'c-');
    %     xlim([min(PPFDGF) max(PPFDGF)]);
    %     ylim([min(GEP) max(GEP)]);
    %     box on; grid on;
    %     title(sprintf('GPP vs PPFD,  b = %3.1f  %3.1f  %3.1f',bGEP));
    %     subplot('position',[0.55 0.28 0.42 0.19]);
    %     plot(mtR,mcR,'o'); fcDateTick(mtR,'Mo',4,1); xlim([min(mtR) max(mtR)]);
    %     title('cR'); ylim([-1 3]);
    %     subplot('position',[0.55 0.04 0.42 0.19]);
    %     plot(mtGEP,mcGEP,'o'); fcDateTick(mtGEP,'Mo',4,1); xlim([min(mtGEP) max(mtGEP)]);
    %     title('cGPP');
    
end;

%	========================================================================
%	========================================================================

