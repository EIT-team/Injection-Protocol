function [w] = create_weight_matrix(vm, Mesh, mat_ref, ROI)
%Calculates weights (total current density in ROI for every possible
%injection pair)
%
% input: vm      potential at all nodes when injection from all electrodes
%                to electrode 1
%        Mesh    Need to have Mesh.sgrad. Can calculate using
%                shape_functions_constructor.m with output of D
%        mat_ref conductivity value in each element 
%        ROI     elements in the ROI
%
% output: w      weight matrix size(n_elec x n_elec) containing total current density in ROI   
%                for every possible injection pair 
%     


%First calculate current density for all electrodes to first electrode
n = 115;
i = 1;
w = zeros(n,n);
[j] = compute_current_density(vm, Mesh, mat_ref, ROI);
j_sum = sum(j);
w(i,i+1:n) = j_sum;

clear j_sum

%Then iterate over all possible combinations
for i = 2:114
    [v] = compute_all_potentials(vm, i, (i-1));
    [j] = compute_current_density(v, Mesh, mat_ref, ROI);
    j_sum = sum(j);
    w(i,i+1:n) = j_sum;
    display(['loop ' num2str(i) ' complete']);
    
    clear v j_sum
end

