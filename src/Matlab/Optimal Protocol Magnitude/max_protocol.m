function [prot] =  max_protocol(w)

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


