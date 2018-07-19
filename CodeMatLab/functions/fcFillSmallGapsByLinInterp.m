	function xGF=fcFillSmallGapsByLinInterp(x,GapSizeMax); 
	
%	xGF=FillSmallGapsByLinInterp(x,GapSizeMax); 
%	
%	FillSmallGapsByLinInterp fills small gaps 
%	in time series by linear interpolation. 
%
%	Gaps are filled if the gap size is less than 
%	or equal to GapSizeMax.

%	Written by Alan Barr 2002.
%	Calls fcGapSize which is also written by Alan Barr 2002.. 

	if sum(~isnan(x))==0; xGF=x; return; end; 
	
	GapSize=fcGapSize(x); 
	iGF=find(GapSize>0 & GapSize<=GapSizeMax); 
	iYaN=find(~isnan(x)); 
	GF=interp1(iYaN,x(iYaN),iGF,'linear'); 
    if isnan(GF)
        GF=0;   % some nan's at start of time series remain, and cause problems downstream
    end
	xGF=x; xGF(iGF)=GF; 	
	