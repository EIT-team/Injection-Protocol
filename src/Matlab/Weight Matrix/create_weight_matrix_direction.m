function [w_mag, w_x, w_y, w_z] = create_weight_matrix_direction(vm, Mesh, mat_ref, ROI)
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


%First calculate current density for 
n = 32;
i = 1;
w_mag = zeros(n,n);
w_x = zeros(n,n);
w_y = zeros(n,n);
w_z = zeros(n,n);

[j_mag, jx, jy, jz] = compute_current_density(vm, Mesh, mat_ref, ROI);
j_sum = sum(j_mag);
j_x = sum(jx);
j_y = sum(jy);
j_z = sum(jz);

w_mag(i,i+1:n) = j_sum;
w_x(i,i+1:n) = j_x;
w_y(i,i+1:n) = j_y;
w_z(i,i+1:n) = j_z;

clear j_sum j_x j_y j_z

for i = 2:n-1
    [v] = compute_all_potentials(vm, i, (i-1), n);
    [j_mag, jx, jy, jz] = compute_current_density(v, Mesh, mat_ref, ROI);
    j_sum = sum(j_mag);
    j_x = sum(jx);
    j_y = sum(jy);
    j_z = sum(jz);
    
    w_mag(i,i+1:n) = j_sum;
    w_x(i,i+1:n) = j_x;
    w_y(i,i+1:n) = j_y;
    w_z(i,i+1:n) = j_z;
    display(['loop ' num2str(i) ' complete']);
    
    clear v j_sum j_x j_y j_z
end

