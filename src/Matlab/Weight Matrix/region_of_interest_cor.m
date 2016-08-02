function [y, Rot] = region_of_interest_cor( Mesh, mat_ref)

test_new = zeros(size(mat_ref,1),1);
x = zeros(size(mat_ref,1),1);
%%To create spherical region of interest%%
%{
ROI_radius = 0.0012;
for i=1:size(mat_ref,1)
    element_centre = mean(Mesh.vtx(Mesh.tri(i,:),:));
    if norm((element_centre-ROI_centre),2) < ROI_radius
        test_new(i) = 1;
    end
end
%}

%%To create a cylindrical region of interest%%

% x = left/right
% y = ventral/dorsal (y in rotated frame goes from -ve to +ve)
% z = posterior/anterior
%ROI_radius = 0.0018;
ROI_radius = 0.00165;
%before z = 0.0183
ROI_centre = [0.02216, 0.00871, 0.01925];
base_ref = [0.0165, 0.01685, 0.01925];
top_ref = [0.0165, 0.00871, 0.01925];

r = ROI_centre - base_ref;
r_ref = top_ref - base_ref;

cos_theta = dot(r, r_ref)/(norm(r,2)*norm(r_ref,2));
theta = acosd(cos_theta);

Rot = [cos_theta, sind(theta); -sind(theta), cos_theta];

ROI_centre_rot = [Rot*ROI_centre(1:2)'; ROI_centre(3)]';

for i=1:size(mat_ref,1)
    element_centre = mean(Mesh.vtx(Mesh.tri(i,:),:));
    element_centre_rot = [Rot*element_centre(1:2)'; element_centre(3)]';
    
    if (element_centre_rot(1) - ROI_centre_rot(1))^2 + (element_centre_rot(3) - ROI_centre_rot(3))^2 < ROI_radius*ROI_radius 
        if element_centre_rot(2) < (ROI_centre_rot(2) + 0.00320)
          if mat_ref(i) == 0.3;
        test_new(i) = 1;
        x(i) = i;
          end
        end
    end
end

y = unique(x);
y = y(2:end);









