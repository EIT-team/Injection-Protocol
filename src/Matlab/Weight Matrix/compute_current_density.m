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


sigma(1:3:3*size(mat_ref,1)) = mat_ref;
sigma(2:3:3*size(mat_ref,1)) = mat_ref;
sigma(3:3:3*size(mat_ref,1)) = mat_ref;

jm=nan(size(ROI,1), size(v,2));
jx=jm;
jy=jm;
jz=jm;

for iInj=1:size(v,2)
    
    Dv = D*v(:,iInj);
    j = Dv.*sigma';
    j2 = [j(1:3:end),j(2:3:end),j(3:3:end)];
    j2 = j2(ROI,:);
    
    jx(:,iInj) = j2(:,1);
    jy(:,iInj) = j2(:,2);
    jz(:,iInj) = j2(:,3);
    
    for j=1:size(j2,1)
        jm(j,iInj) = norm(j2(j,:),2);
    end
    
    clear Dv j j2
    
end
toc
end

