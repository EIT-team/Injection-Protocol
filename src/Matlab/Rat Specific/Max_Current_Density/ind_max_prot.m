function [prot] = ind_max_prot(prot_max, n)



prot=[];
elec = zeros(114,1); 
while size(prot_max,1) > 0
    prot = [prot; prot_max(1,:)];
    
    prot_max(1,:) = 0;
    
    for i = 1:size(prot,1)
        a = prot(i,1);
        b = prot(i,2);
        
        elec(a) = elec(a) + 1;
        elec(b) = elec(b) + 1;
    end
    
    A = find(elec >= n);
    
    for j = 1:size(A,1)
        for i = 1:size(prot_max,1)
            if prot_max(i,1) == A(j) || prot_max(i,2) == A(j)
             prot_max(i,:) = 0;
            end
        end
    
    end
     
   B = find(prot_max(:,1) ~=0);
   prot_max = prot_max(B,:);
   elec = zeros(114,1);
    
end

    
        
       