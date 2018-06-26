function [ROI, ROI_nodes] = region_of_interest_sphere(Mesh,ROI_centre,ROI_radius)
%Finds elements that are in a specified region of interest in brain
%Finds elements that lie in a sphere about ROI centre

%Input
%Mesh           mesh that you are using

%Output
%ROI            elements contained within spherical region of interest
%ROI nodes      nodes of all elements (easier for visualisation check in matlab)

% x = left/right
% y = ventral/dorsal (y in rotated frame goes from -ve to +ve)
% z = posterior/anterior

%Centre for VPM
%ROI_centre = [0.0192, 0.01430, 0.01801];

%Centre and radius for VPL and PO
%Right hemisphere
%ROI_centre = [0.0193, 0.0131, 0.01801];
%Left hemisphere
%ROI_centre = [0.01392, 0.0131, 0.01801];
%ROI_radius = 0.0015;

%Centre and radius for VPL
%Run program twice with the two different ROI_centre and take setdiff to 
%get elements in a crescent
%Right hemisphere
%ROI_centre = [0.0193, 0.0133, 0.01811];
%ROI_centre = [0.0184,0.0129,0.01781];
%Left hemisphere
%ROI_centre = [0.01392, 0.0133, 0.01811];
% ROI_centre = [0.01482, 0.0129, 0.01781];
% ROI_radius = 0.0015;

%For cylindrical mesh

% ROI_centre = [150,100,150]/1000;
% ROI_radius = 10/1000;

ROI = zeros(size(Mesh.tri,1),1);


cnts=(Mesh.vtx(Mesh.tri(:,1),:)+Mesh.vtx(Mesh.tri(:,2),:)+Mesh.vtx(Mesh.tri(:,3),:)+Mesh.vtx(Mesh.tri(:,4),:))./4;



for iElem=1:size(Mesh.tri,1)
    element_centre = cnts(iElem,:);
    if norm((element_centre-ROI_centre),2) < ROI_radius
        ROI(iElem) = iElem;
    end
end

ROI = unique(ROI);
ROI = ROI(2:end);

ROI_nodes = Mesh.tri(ROI);
ROI_nodes = reshape(ROI_nodes,size(ROI_nodes,1)*size(ROI_nodes,2),1);
ROI_nodes = unique(ROI_nodes);

end





