function [w_mag, w_x, w_y, w_z] = create_weight_matrix_direction(vm, Mesh, mat_ref, ROI,n_elec)
%Calculates weights (total current density in ROI for every possible
%injection pair)
%
% input: vm      potential at all nodes when injection from all electrodes
%                to electrode 1 (get this from Fwd.current_field(end-nodenum+1:end,:);
%        Mesh    Need to have Mesh.sgrad. Can calculate using
%                shape_functions_constructor.m with output of D
%        mat_ref conductivity value in each element 
%        ROI     elements in the ROI
%
% output: w      weight matrix size(n_elec x n_elec) containing total current density in ROI   
%                for every possible injection pair 
%     


%First calculate current density for 
% n_elec = 32;
iInj = 1;
w_mag = zeros(n_elec,n_elec);
w_x = zeros(n_elec,n_elec);
w_y = zeros(n_elec,n_elec);
w_z = zeros(n_elec,n_elec);

[j_mag, jx, jy, jz] = compute_current_density(vm, Mesh, mat_ref, ROI);
j_sum = sum(j_mag);
j_x = sum(jx);
j_y = sum(jy);
j_z = sum(jz);

w_mag(iInj,iInj+1:n_elec) = j_sum;
w_x(iInj,iInj+1:n_elec) = j_x;
w_y(iInj,iInj+1:n_elec) = j_y;
w_z(iInj,iInj+1:n_elec) = j_z;

clear j_sum j_x j_y j_z

for iInj = 2:n_elec-1
    [v] = compute_all_potentials(vm, iInj, (iInj-1), n_elec);
    [j_mag, jx, jy, jz] = compute_current_density(v, Mesh, mat_ref, ROI);
    j_sum = sum(j_mag);
    j_x = sum(jx);
    j_y = sum(jy);
    j_z = sum(jz);
    
    w_mag(iInj,iInj+1:n_elec) = j_sum;
    w_x(iInj,iInj+1:n_elec) = j_x;
    w_y(iInj,iInj+1:n_elec) = j_y;
    w_z(iInj,iInj+1:n_elec) = j_z;
    disp(['loop ' num2str(iInj) ' complete']);
    
    clear v j_sum j_x j_y j_z
end

