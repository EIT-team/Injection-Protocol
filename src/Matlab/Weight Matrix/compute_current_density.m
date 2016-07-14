function [jm, jx, jy, jz] = compute_current_density(v, Mesh, mat_ref, ROI)
% Computes the current density on a mesh based on the potential
% distribution
%
% input: mesh       mesh with vtx and tri and mat_ref (struct)
%        v          potentials on all nodes of the mesh (vector)
%        mat_ref    conductivity distribution (vector)
%        ROI        element indices of ROI

% output: jm        current density vector for each element in ROI
%         jx        current density vector in x direction for each element
%         jy        current density vector in y direction for each element
%         jz        current density vector in z direction for each element 
%  


tic
D = Mesh.sgrad;
Dv = D*v;

sigma(1:3:3*size(mat_ref,1)) = mat_ref;
sigma(2:3:3*size(mat_ref,1)) = mat_ref;
sigma(3:3:3*size(mat_ref,1)) = mat_ref;


for i=1:size(v,2)
j(:,i) = Dv(:,i).*sigma';
end

for i=1:size(v,2)
j2{i} = [j(1:3:end,i),j(2:3:end,i),j(3:3:end,i)];
end

for i = 1:size(v,2)
    j2{i} = j2{i}(ROI,:);
end

for i = 1:size(v,2)
    jx(:,i) = j2{i}(:,1);
    jy(:,i) = j2{i}(:,2);
    jz(:,i) = j2{i}(:,3);
end

for k=1:size(v,2)
    for i=1:size(j2{1,1},1)
        jm(i,k) = norm(j2{1,k}(i,:),2);
    end
end
toc

