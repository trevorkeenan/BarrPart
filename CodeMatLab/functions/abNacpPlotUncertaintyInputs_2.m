	function abNacpPlotUncertaintyInputs_2(t,NEE,uStar,PPFD,Ta,Ts,Vpd,fNight,cSiteYr,iFig,cFlux) 

% abNacpPlotUncertaintyInputs
%	plots annual time series for input data 
%	in NACP NEE uncertainty analysis.

    if length(t)<=8784
        nRecsPerDay=24;
    else
        nRecsPerDay=48;
    end
% 	nRecsPerDay=round(1/nanmedian(diff(t)));
  	ntNight=sum(reshape(fNight,nRecsPerDay,length(fNight)/nRecsPerDay));
	mtNight=mean(reshape(t,nRecsPerDay,length(t)/nRecsPerDay));
	
	fcFigLoc(iFig,0.6,0.9,'SE');
	
	for iPlot=1:7; 
		
		switch iPlot; 
			
			case 1; p=[0.05 0.70 0.42 0.25]; x=t; y=[Ta Ts]; cy=[cSiteYr '  Ta Ts']; 
			case 2; p=[0.05 0.375 0.42 0.25]; x=t; y=Vpd; cy='VPD'; 
			case 3; p=[0.05 0.05 0.42 0.25]; x=t; y=PPFD; cy='PPFD'; 
			
			case 4; p=[0.55 0.70 0.42 0.25]; x=t; y=uStar; cy='uStar'; 
			case 5; p=[0.55 0.375 0.42 0.25]; x=t; y=NEE; cy=cFlux;
			case 6; p=[0.55 0.05 0.42 0.25]; x=mtNight; y=ntNight; cy='nPeriods Per Night';
			
		end; 
		
		subplot('position',p);
		plot(x,y,'.'); fcDateTick(t,'Mo',4,1); xlim([min(t) max(t)]);
		title(sprintf('%s   nMiss %g %g',cy,sum(isnan(y))));
		
	end; 		
			
	
