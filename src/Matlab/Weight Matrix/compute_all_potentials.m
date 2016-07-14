function [v] = compute_all_potentials(vm, n, i, el)
%Calculates potential resulting from all possible current injection pairs
%
% input: vm      potential at all nodes when injection from all electrodes
%                to electrode 1
%        n       starting value of iteration
%        i       the electrode you want to calculate potential with respect
%                to e.g electrode 57
%        el      total number of electrodes
%
% output: v      potential resulting from all combinations of electrodes
%                between reference and others 
%     

for k = n:el-1;
v(:,k-(n-1)) = vm(:,i) - vm(:,k); %e.g 2 -> 3, 2 -> 4
end

