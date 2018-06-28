function [Protocol] = create_protocol(w,n)
%Gives out protocol that maximises the weight while finding independent
%injection pairs
%
% input:    w       matrix of weights (input the inverse to find maximum)
% output:   protocol    optimal protocol
% pos : electrode positioning nx3

%    w = w(1:114,1:114);
%    for i = 1:114
%        w(i:114,i) = w(i,i:114);
%    end
   
    w = w(1:n,1:n);
   for i = 1:n
       w(i:n,i) = w(i,i:n);
   end
   
   w = 1./w;

   T=tril(sparse(w));
   [ST] = graphminspantree(T);
   [i,j,s]=find(ST);
   [A]=sortrows([s,i,j]);
   Protocol = A(:,2:3);
