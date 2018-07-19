% this script will...
%   loop through each site
%   and site year
%   and pass the current data to Alan Barrs uncertainty script

%               RUNs IN MATLAB 2017a

% this is the same as uncertaintyMaster.m
% except the uStar threshold is set to PI preferred values

% processing time: ~30s per year for nBoot = 2

% written by T. Keenan as a wrapper to the BarrPart gap-filling/partitioning code of Alan Barr


% addpath(genpath('./functions/'),'./functions2/');
addpath(genpath('./functions/'));
siteResultsLoc = '../../../dataPartitioned/';

% set the number of bootstrap samples
%nBoot=2; % nBoot=1 may cause problems
nBoot=200;

% 1. load the uStar thresholds for each site
% uStars=importdata('../../uStarThresholds.csv');
ustarthreshold=0.5; %uStars.data;

% set data column locations
columns.year=1;
columns.ddoy=2;

columns.NEE=5;
columns.LE=7;
columns.H=6;

columns.ustar=8;
columns.Tair=3;
columns.SoilT=3;
columns.PAR=4;
columns.VPD=3; % NOTE: VPD NOT USED 

% columns of the gap filled output and uncertainties
columns.NEEf=9;
columns.NEEu=10;
columns.GPPf=11;
columns.GPPu=12;
columns.REf=13;
columns.REu=14;
columns.LEf=15;
columns.LEu=16;
columns.Hf=17;
columns.Hu=18;

% set data column of the output data
% these are same location as above
% but can have different names
columns.out.year=1;
columns.out.ddoy=2;
columns.out.TA=3;
columns.out.PPFDin=4;
columns.out.FC=5;
columns.out.H=6;
columns.out.LE=7;
columns.out.ustar=8;
columns.out.NEEf=9;
columns.out.NEEu=10;
columns.out.GPPf=11;
columns.out.GPPu=12;
columns.out.REf=13;
columns.out.REu=14;
columns.out.LEf=15;
columns.out.LEu=16;
columns.out.Hf=17;
columns.out.Hu=18;

% % gap filled met columns
% columns.PREC_f=19;
% columns.Rgl_f=20;
% columns.PRESS_f=21;
% columns.WS_f=22;
% columns.SWC1_f=23;
% columns.Rg_f=24;
% columns.RH_f=25;


% loop throught each site
% extract data for each year
% and send it for BARR uncertainty analysis
siteFoldersLoc = '../../../dataToPartition/';
files=dir(siteFoldersLoc);
% remove '.','..','./DS_store'
files = files(~strncmpi('.', {files.name}, 1));

% initialize the results matrix
randomU=cell(length(files),20);  % not sure what the 20 is here
ustarU=cell(length(files),20);

%for i=1:length(ustarthreshold)
for i=1:1
    clc
    
    % set the current site name
    cSite=files(i).name;
    
    % identify the number of years in the current site
    years=dir(strcat(siteFoldersLoc,cSite,'/*.csv'));
    
    numYears=length(years);
    
    disp(cSite);
    
    for j= 1:numYears
        %     for j=1:numYears
        
        tic;
        % Open the text file.
        filename=strcat(siteFoldersLoc,cSite,'/',years(j).name);
        dataArray = dlmread(filename,',',1,0);
        
        % if there is no uStar info, there should be no flux data
        indX=dataArray(:,columns.ustar)<=-999;
        dataArray(indX,columns.LE)=NaN;
        dataArray(indX,columns.H)=NaN;
        dataArray(indX,columns.NEE)=NaN;
        dataArray(indX,columns.ustar)=NaN;
        
        dataL2_4alan=dataArray;
        
        % get the current year
        cYr=num2str(dataArray(2,1));
        disp(cYr);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        [randomU_LE{i,j}, ustarU_LE{i,j},LEQc,LEgf,bNegLE,bPosLE] ...
            = abScriptNacpUncertaintyFcrn_fixedUstarDP_LE(dataL2_4alan,cSite,cYr,columns,nBoot,ustarthreshold(i));
        
        % gap fill the H fluxes
        [randomU_H{i,j}, ustarU_H{i,j},HQc,Hgf,bNegH,bPosH] ...
            = abScriptNacpUncertaintyFcrn_fixedUstarDP_H(dataL2_4alan,cSite,cYr,columns,nBoot,ustarthreshold(i));
        
        % gap fill and partition the carbon fluxes
        [randomU{i,j}, ustarU{i,j},NEEQc,NEEgfgC,bNegNEE,bPosNEE] ...
            = abScriptNacpUncertaintyFcrn_fixedUstarDP(dataL2_4alan,cSite,cYr,columns,nBoot,ustarthreshold(i));
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fill the data array with the gap-filled data
        
        dataArray2=dataArray;
        
        % NEE
        dataArray2(:,columns.NEEf)=-1*NEEgfgC;
        
        tmp=zeros(1,length(NEEQc))+NaN;
        % for negative data the error follows bNeg
        indX3=NEEgfgC<0;
        tmp(indX3)=bNegNEE(1)+bNegNEE(2)*NEEgfgC(indX3);
        
        % for positive original data the error follows bPos
        indX3=NEEgfgC>=0;
        tmp(indX3)=bPosNEE(1)+bPosNEE(2)*NEEgfgC(indX3);
        
        dataArray2(:,columns.NEEu)=tmp;
        
        % GPP
        dataArray2(:,columns.GPPf)=mean(randomU{i,j}.GPP.Hourly);
        dataArray2(:,columns.GPPu)=std(randomU{i,j}.GPP.Hourly);
        
        % RE
        dataArray2(:,columns.REf)=mean(randomU{i,j}.RE.Hourly);
        dataArray2(:,columns.REu)=std(randomU{i,j}.RE.Hourly);
        
        % LE
        dataArray2(:,columns.LEf)=LEgf;
        
        % estimate the error
        tmp=bPosLE(1)+bPosLE(2)*LEgf;
        dataArray2(:,columns.LEu)=tmp;
        
        % H
        dataArray2(:,columns.Hf)=Hgf;
        
        % estimate the error
        tmp=zeros(1,length(Hgf))+NaN;
        % for negative data the error follows bNeg
        indX3=Hgf<0;
        tmp(indX3)=bPosH(1)+bPosH(2)*Hgf(indX3);
        
        % for positive original data the error follows bPos
        indX3=Hgf>=0;
        tmp(indX3)=bPosH(1)+bPosH(2)*Hgf(indX3);
        
        dataArray2(:,columns.Hu)=tmp;
        
        dataArray2(isnan(dataArray2))=-9999;
        dataArray2(isinf(dataArray2))=-9999;
        
        % write the results to file
        names = fieldnames(columns.out);
        outTable = array2table(dataArray2,'VariableNames',names);
        
        filename = strcat(cSite,'/HalfHourly/',cSite,'_',num2str(cYr),'.csv');
        writetable(outTable,strcat(siteResultsLoc,filename))
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         %%          DAILY
        %
        %         % 1. Aggregate met data to Daily (currently T, VPD, PAR)
        dataArray2(dataArray2==-9999)=NaN;
        %
        %         % identify the days
        allDays = floor(dataArray2(:,2));
        days=unique(allDays);
        %
        %         % loop through the days
        daily=zeros(365,4)+NaN;
        for ii=1:length(days)
            currentDOY=days(ii);
            
            % find data location for current DOY
            indX=allDays==currentDOY;
            daily(ii,1)=mean(dataArray2(indX,1)); % Year
            daily(ii,2)=currentDOY;               % DOY
            daily(ii,3)=mean(dataArray2(indX,columns.Tair)); % Air Temperature
            %             daily(ii,4)=mean(dataArray2(indX,columns.SoilT)); % Soil Temperature
            daily(ii,4)=sum(dataArray2(indX,columns.PAR)); % PAR
            %             daily(ii,6)=mean(dataArray2(indX,columns.VPD)); % VPD
            
            %             daily(ii,7)=nansum(dataArray2(indX,columns.PREC_f)); % Precip
            %             daily(ii,8)=nansum(dataArray2(indX,columns.Rgl_f)); %
            %             daily(ii,9)=nanmean(dataArray2(indX,columns.PRESS_f)); % Pressure
            %             daily(ii,10)=nanmean(dataArray2(indX,columns.WS_f));
            %             daily(ii,11)=nanmean(dataArray2(indX,columns.SWC1_f));
            %             daily(ii,12)=nansum(dataArray2(indX,columns.Rg_f));
            %             daily(ii,13)=nanmean(dataArray2(indX,columns.RH_f));
            
        end
        
        % get the flux data and uncertainties
        clear GF_fixedUstar
        tmp=length(nanmean(randomU{i,j}.NEE.Daily));
        GF_fixedUstar(1,1:tmp)=str2double(cYr);
        GF_fixedUstar(2,1:tmp)=mean(randomU{i,j}.NEE.Daily);
        GF_fixedUstar(3,1:tmp)=std(randomU{i,j}.NEE.Daily);
        GF_fixedUstar(4,1:tmp)=mean(randomU{i,j}.GPP.Daily);
        GF_fixedUstar(5,1:tmp)=std(randomU{i,j}.GPP.Daily);
        GF_fixedUstar(6,1:tmp)=mean(randomU{i,j}.RE.Daily);
        GF_fixedUstar(7,1:tmp)=std(randomU{i,j}.RE.Daily);
        GF_fixedUstar(8,1:tmp)=mean(randomU_LE{i,j}.NEE.Daily);
        GF_fixedUstar(9,1:tmp)=std(randomU_LE{i,j}.NEE.Daily);
        GF_fixedUstar(10,1:tmp)=mean(randomU_H{i,j}.NEE.Daily);
        GF_fixedUstar(11,1:tmp)=std(randomU_H{i,j}.NEE.Daily);
        
        tmp=horzcat(daily,GF_fixedUstar(2:end,:)');
        tmp(isnan(tmp))=-9999;
        
        
        % write the results to file
        names = fieldnames(columns.out);
        names2 = names([1,2,3,4,9:end],1);
        outTable = array2table(tmp,'VariableNames',names2);
        
        filename = strcat(cSite,'/Daily/',cSite,'_',num2str(cYr),'.csv');
        writetable(outTable,strcat(siteResultsLoc,filename))
        
        
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         %%               Monthly
        %%
        % 1. Aggregate met data to Monthly (currently T, VPD, PAR)
        
        % loop through the days
        monthly=zeros(12,4)+NaN;
        monthx=1;
        for ii=1:30:(length(daily)-30)
            
            monthly(monthx,1)=mean(daily(ii:ii+29,1)); % Year
            monthly(monthx,2)=monthx;               % Month
            monthly(monthx,3)=mean(daily(ii:ii+29,3)); % Air Temperature
            %             monthly(monthx,4)=mean(daily(ii:ii+29,4)); % Soil T
            monthly(monthx,4)=sum(daily(ii:ii+29,4)); %
            %             monthly(monthx,6)=mean(daily(ii:ii+29,6)); % VPD
            
            %             monthly(monthx,7)=sum(daily(ii:ii+29,7));
            %             monthly(monthx,8)=sum(daily(ii:ii+29,8));
            %             monthly(monthx,9)=mean(daily(ii:ii+29,9));
            %             monthly(monthx,10)=mean(daily(ii:ii+29,10));
            %             monthly(monthx,11)=mean(daily(ii:ii+29,11));
            %             monthly(monthx,12)=sum(daily(ii:ii+29,12));
            %             monthly(monthx,13)=mean(daily(ii:ii+29,13));
            monthx=monthx+1;
        end
        
        clear GF_fixedUstar
        GF_fixedUstar(1,1:12)=str2double(cYr);
        GF_fixedUstar(2,1:12)=mean(randomU{i,j}.NEE.Monthly);
        GF_fixedUstar(3,1:12)=std(randomU{i,j}.NEE.Monthly);
        GF_fixedUstar(4,1:12)=mean(randomU{i,j}.GPP.Monthly);
        GF_fixedUstar(5,1:12)=std(randomU{i,j}.GPP.Monthly);
        GF_fixedUstar(6,1:12)=mean(randomU{i,j}.RE.Monthly);
        GF_fixedUstar(7,1:12)=std(randomU{i,j}.RE.Monthly);
        GF_fixedUstar(8,1:12)=mean(randomU_LE{i,j}.NEE.Monthly);
        GF_fixedUstar(9,1:12)=std(randomU_LE{i,j}.NEE.Monthly);
        GF_fixedUstar(10,1:12)=mean(randomU_H{i,j}.NEE.Monthly);
        GF_fixedUstar(11,1:12)=std(randomU_H{i,j}.NEE.Monthly);
        
        
        tmp=horzcat(monthly,GF_fixedUstar(2:end,:)');
        tmp(isnan(tmp))=-9999;
        outTable = array2table(tmp,'VariableNames',names2);
        
        filename = strcat(cSite,'/Monthly/',cSite,'_',num2str(cYr),'.csv');
        writetable(outTable,strcat(siteResultsLoc,filename))
        %%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% ANNUAL
        
        clear GF_fixedUstar
        fname=strcat('../../../dataPartitioned/',cSite,'/Annual/allYears_f.mat');
        
        if exist(fname,'file')
            load(fname)
        end
        
        GF_fixedUstar(j,1)=str2double(cYr);
        %         indX=dataArray(:,1)==GF_fixedUstar(j,1);
        
        GF_fixedUstar(j,2)=mean(dataArray2(:,columns.Tair)); % Air T
        
        %         GF_fixedUstar(j,3)=mean(dataArray2(:,columns.SoilT)); % Soil Temperature
        GF_fixedUstar(j,3)=sum(dataArray2(:,columns.PAR)); % PAR
        %         GF_fixedUstar(j,5)=mean(dataArray2(:,columns.VPD)); % VPD
        
        %         GF_fixedUstar(j,6)=sum(dataArray2(:,columns.PREC_f)); % Precip
        %         GF_fixedUstar(j,7)=sum(dataArray2(:,columns.Rgl_f)); %
        %         GF_fixedUstar(j,8)=mean(dataArray2(:,columns.PRESS_f)); % Pressure
        %         GF_fixedUstar(j,9)=mean(dataArray2(:,columns.WS_f));
        %         GF_fixedUstar(j,10)=mean(dataArray2(:,columns.SWC1_f));
        %         GF_fixedUstar(j,11)=sum(dataArray2(:,columns.Rg_f));
        %         GF_fixedUstar(j,12)=mean(dataArray2(:,columns.RH_f));
        
        GF_fixedUstar(j,4)=mean(randomU{i,j}.NEE.Annual);
        GF_fixedUstar(j,5)=std(randomU{i,j}.NEE.Annual);
        GF_fixedUstar(j,6)=mean(randomU{i,j}.GPP.Annual);
        GF_fixedUstar(j,7)=std(randomU{i,j}.GPP.Annual);
        GF_fixedUstar(j,8)=mean(randomU{i,j}.RE.Annual);
        GF_fixedUstar(j,9)=std(randomU{i,j}.RE.Annual);
        GF_fixedUstar(j,10)=mean(randomU_LE{i,j}.NEE.Annual);
        GF_fixedUstar(j,11)=std(randomU_LE{i,j}.NEE.Annual);
        GF_fixedUstar(j,12)=mean(randomU_H{i,j}.NEE.Annual);
        GF_fixedUstar(j,13)=std(randomU_H{i,j}.NEE.Annual);
        
        
        save(fname, 'GF_fixedUstar')
        
        elapsedTime=toc;
        disp('Processing time per year:')
        disp(elapsedTime)
        
    end
    
    % save the annual csv file after processing all years
    clear GF_fixedUstar
    fname=strcat('../../../dataPartitioned/',cSite,'/Annual/allYears_f.mat');
    
    if exist(fname,'file')
        load(fname)
    end
    
    tmp=GF_fixedUstar;
    filename=strcat(cSite,'/Annual/',cSite,'allYears_af.csv');
    
    names3 = names2([1,3:end],1);
     
    outTable = array2table(tmp,'VariableNames',names3);
    
    writetable(outTable,strcat(siteResultsLoc,filename))
    
    delete(fname)
    
end


