function [prot,cd] =  sort_protocol(w)
%sort_protocol - sort all possible protocol lines in order of maximum current density in ROI

%Input w        weight matrix with current density in ROI for every
%               possible injection pair
%output prot    injection pairs in order of max current density
%       cd      corresponding current density for each pair



iInj=1;
T = triu(w);

nInj=nnz(w);
prot=zeros(nInj,2);
cd=zeros(nInj,1);
while sum(T(:)) ~= 0
    [maxcd, loc] = max(T(:));
    [R,C] = ind2sub(size(T),loc);
    prot(iInj,1) = R;
    prot(iInj,2) = C;
    cd(iInj)=maxcd;
    
    T(R,C) = 0;
    iInj = iInj+ 1;
end
end