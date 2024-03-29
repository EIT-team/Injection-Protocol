function [maxprot,prot_all,cdd] = max_parallel(w,nInj)
%MAX_PARALLEL Finds protocol with highest current density in ROI with no
%repeated electrodes - for parallel injection
%%
[protfull,cd] =  sort_protocol(w);
%% find the current density for all possible pairs

nRows=size(protfull,1);
cdd=zeros(nRows,1);
prot_idx=zeros(nRows,nInj);
prot_all=zeros(nInj,2,nRows);

% for each potential injection, find the best non repeating pairs
for iRow=1:nRows
    
    pout = findInjs(iRow,protfull,nInj);
    pout=sortrows(pout);
    
    prot_idx(iRow,:)=find(ismember(protfull,pout,'rows'));
    
    cdd(iRow)=sum(cd(prot_idx(iRow,:)));
    
    
    prot_all(:,:,iRow)=pout;
    
end
%% find the best combination

[maxcd,maxidx]=max(cdd);

maxprot=sortrows(protfull(prot_idx(maxidx,:)',:));

%% remove non unique ones

%sort first
[cdd,c_idx]=sort(cdd,'descend');
prot_all=prot_all(:,:,c_idx);


[cdd, c_idx]=(unique(cdd,'stable'));
prot_all=prot_all(:,:,c_idx);


end

function pout = findInjs(startpos,protfull,nInj)
%finds the set of injection pairs which gives the highest current density
%given a starting row

prot=protfull; % potential injections
pout=zeros(nInj,2);
%first one is always the starting point
pout(1,:)=prot(startpos,:);
prot(startpos,:)=[];


for iInj = 2:nInj
    %find protocol lines which have only unusued electrodes
    prot_cand=~any(ismember(prot,pout),2);
    %as list is sorted, we can take the first unused one as this will be
    %the max
    Inj_idx=find(prot_cand,1);
    pout(iInj,:)=prot(Inj_idx,:);
end

end
