	function [gs]=fcGapSize(x); 
	
%	[gs]=fcGapSize(x); 
%	
%	gapsize determines the number of contiguous missing data 
%	for each element of the column vector x.
	
%	Written by Alan Barr 2002.
   
	fNaN=isnan(x); gs=zeros(size(x)); lx=length(x); 
	nNaN=sum(fNaN); if nNaN==0; return; end; 
	iGapStart=find(diff(fNaN)==1)+1; iGapStart=fcx2colvec(iGapStart); 
		if fNaN(1)==1; iGapStart=[1; iGapStart]; end; 
	iGapEnd=find(diff(fNaN)==-1); iGapEnd=fcx2colvec(iGapEnd); 
		if fNaN(lx)==1; iGapEnd=[iGapEnd; lx]; end; 
	nGaps=length(iGapStart); 
	for i=1:nGaps; 
		gs(iGapStart(i):iGapEnd(i))=iGapEnd(i)-iGapStart(i)+1;
	end;    
