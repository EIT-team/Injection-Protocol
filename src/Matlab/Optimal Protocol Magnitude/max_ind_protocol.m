function [Protocol] =  max_ind_protocol(w)
%%Finds independent protocol by that maximises current density in ROI but
%%only uses each electrode twice

%Input w        weight matrix with current density in ROI for every
%               possible injection pair
%
%Output Protocol protocol that uses each electrode twice


%%
%First sort the protocol lines in order of maximum current density of ROI
i=1;
T = triu(w);
while sum(T(:)) ~= 0
    [~, loc] = max(T(:));
    [R,C] = ind2sub(size(T),loc);
    prot(i,1) = R;
    prot(i,2) = C;
    
    T(R,C) = 0;
    i = i+ 1;
end


%%
%Then go through each one and select the protocol line that maximises
%current density but doesnt use an electrode more than twice
n = size(w,1);

Protocol=[];
elec = zeros(n,1);

while size(prot,1) > 0
    Protocol = [Protocol; prot(1,:)];
    
    prot(1,:) = 0;
    
    for i = 1:size(Protocol,1)
        a = Protocol(i,1);
        b = Protocol(i,2);
        
        elec(a) = elec(a) + 1;
        elec(b) = elec(b) + 1;
    end
    
    A = find(elec >= 2);
    
    for j = 1:size(A,1)
        for i = 1:size(prot,1)
            if prot(i,1) == A(j) || prot(i,2) == A(j)
                prot(i,:) = 0;
            end
        end
        
    end
    
    B = find(prot(:,1) ~=0);
    prot = prot(B,:);
    elec = zeros(n,1);
    
end
end