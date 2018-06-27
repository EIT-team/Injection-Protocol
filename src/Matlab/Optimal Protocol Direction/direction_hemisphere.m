<<<<<<< HEAD:src/Matlab/Optimal Protocol Direction/angle_with_order.m
%%
clear all
load('..\..\..\resources\data\Direction weight matrix\direction_rat_VPL_LH.mat');

start_method = 'A'; %A for first start search, B for alternative start of search
%A good way to decide which option is better is to look at total_j and
%zero_bins - want total_j to be as big as possible but don't want too many
%zero bins

%%
% n = 32;

% M = w_mag(1:n,1:n);
% X = w_x(1:n,1:n);
% Y = w_y(1:n,1:n);
% Z = w_z(1:n,1:n);

n = 114;
=======
n = 52;
>>>>>>> master:src/Matlab/Optimal Protocol Direction/direction_hemisphere.m

% M = EX.w_mag(1:n,1:n);
% X = EX.w_x(1:n,1:n);
% Y = EX.w_y(1:n,1:n);
% Z = EX.w_z(1:n,1:n);

M = w_mag(1:n,1:n);
X = w_x(1:n,1:n);
Y = w_y(1:n,1:n);
Z = w_z(1:n,1:n);

<<<<<<< HEAD:src/Matlab/Optimal Protocol Direction/angle_with_order.m
w_cart=cell(n,n);
rho=zeros(n,n);
phi=zeros(n,n);
theta=zeros(n,n);
w_spher=cell(n,n);

%Convert from cartesian to spherical coordinates
%Take abs of all values as we only want to look in 1/8 of sphere
=======

>>>>>>> master:src/Matlab/Optimal Protocol Direction/direction_hemisphere.m
for j = 1:n
    for k = 1:n
        w_cart{j,k} = [X(j,k),Y(j,k),Z(j,k)];
        rho(j,k) = norm([X(j,k), Y(j,k), Z(j,k)],2);
        
        if X(j,k) < 0
            X(j,k) = -1*X(j,k);
            Y(j,k) = -1*Y(j,k);
            Z(j,k) = -1*Z(j,k);
        phi(j,k) = atan((X(j,k))/(Y(j,k)))*(180/pi);
        theta(j,k) = acos((Z(j,k))/rho(j,k))*(180/pi);
     
        else
        phi(j,k) = atan((X(j,k))/(Y(j,k)))*(180/pi);
        theta(j,k) = acos((Z(j,k))/rho(j,k))*(180/pi);
    end
        
        w_spher{j,k} = [rho(j,k),phi(j,k),theta(j,k)];
    end
end

<<<<<<< HEAD:src/Matlab/Optimal Protocol Direction/angle_with_order.m
%%
%Define the number of angles you want to segment into
%Number of bins will be n_angle*n_angle as we segment in both theta and phi
%direction
n_angles = 9;
d_angle = 90/n_angles;
idx=cell(n_angles,n_angles);
%Find all protocols that lie in each bin
for iphi = 1:n_angles
    for itheta = 1:n_angles
        idx{iphi,itheta} = find(phi >= (iphi-1)*d_angle & phi < iphi*d_angle & theta >= (itheta - 1)*d_angle & theta < itheta*d_angle);
=======

n_angles = 6;
d_angle = 180/n_angles;
%phi goes from 0 to 90 and then -90 to 0
%theta goes from  0 to 180

phi_seg = [0:d_angle:90, -90:d_angle:0];
theta_seg = [0:d_angle:180];

for iphi = 1:length(phi_seg)-1
    for itheta = 1:length(theta_seg)-1
        idx{iphi,itheta} = find(phi >= phi_seg(iphi) & phi < phi_seg(iphi+1) & theta >= theta_seg(itheta) & theta < theta_seg(itheta+1));
>>>>>>> master:src/Matlab/Optimal Protocol Direction/direction_hemisphere.m
    end
end

idx( all(cellfun(@isempty,idx),2), : ) = [];

%%
idx = reshape(idx, n_angles*n_angles,1);
rm = find(cellfun(@isempty,idx));
idx(rm) = [];

n_bins = size(idx,1);

<<<<<<< HEAD:src/Matlab/Optimal Protocol Direction/angle_with_order.m
%%
if start_method == 'A'
    %Find the protocol with the maximum magnitude in each bin and sort from
    %largest to smallest
    %Determines the order with which we look through the bins as later ones
    %will be more constrained as fewer allowable electrodes left
    max_j=zeros(n_bins);
    for iMag = 1:n_bins
        max_j(iMag) = max(M(idx{iMag,1}));
    end
    
    [~,I] = sort(max_j, 'descend');
    
end
%%
if start_method == 'B'
    %This is an alternative option to start the search
    %Finds the total number of possible electrodes that can be used to address
    %a bin
    %The lower the number the fewer possibilities so better to start with these
    %first to ensure there is at least one protocol that targets that bin
    
    elec =cell(n_bins,1);
    
    for i = 1:n_bins
        [R,C] = ind2sub(size(theta), idx{i,1}(:));
        diff_elec = unique([R;C]);
        elec{i,1} = diff_elec;
    end
    
    [~,I] = sort(cellfun(@length,elec));
    
end

%%
%First go through each bin and find the best protocol line, making sure
%they are all independent

%Store value of total current density in each bin
=======

for i = 1:n_bins
    [R,C] = ind2sub(size(theta), idx{i,1}(:));
    diff_elec = unique([R;C]);
    elec{i,1} = diff_elec;
end

[~,I] = sort(cellfun(@length,elec));

>>>>>>> master:src/Matlab/Optimal Protocol Direction/direction_hemisphere.m
j_bins = zeros(n_bins,1);
% prot=zeros(n_bins,2);
for iProt = 1:n_bins
    [max_val,loc] = max(M(idx{I(iProt),1}));
    loca = idx{I(iProt),1}(loc);
    
    %Remove the index we have chosen so we don't consider it next time
    idx{I(iProt),1}(loc) = [];
    [R,C] = ind2sub(size(theta), loca);
    prot(iProt,1) = R;
    prot(iProt,2) = C;
    
    j_bins(I(iProt)) = max_val;
    
    %Keep track of which electrodes have been used and previous protocol lines
    %so we keep protocol independent
    count = zeros(n,1);
    for iCount = 1:size(prot,1)
        a = prot(iCount,1);
        b = prot(iCount,2);
        
        count(a) = count(a) + 1;
        count(b) = count(b) + 1;
        
        M(prot(iCount,1),prot(iCount,2)) = 0;
        M(prot(iCount,2),prot(iCount,1)) = 0;
    end
    
    A = find(count >= 2);
    
    M(A,:) = 0;
    M(:,A) = 0;
end

%%
for iProt = n_bins:n-1
    %    for iProt = 1:n-1
    
    %sort current density from smallest to largest to see which bin we should
    %be adding to
    [~,NI] = sort(j_bins);
    
    %Sometimes possible that we can't add to that smallest bin any more as no indpendent
    %pairs left
    %So loop until what we add is a non-zero value
    for ibin = 1:length(NI)
        [max_val, loc] = max(M(idx{NI(ibin),1}));
        if max_val ~= 0
            %As soon as we find a non-zero value break out from the loop
            n_idx = NI(ibin);
            break;
        end
    end
    
    loca = idx{n_idx,1}(loc);
    idx{n_idx,1}(loc) = [];
    [R,C] = ind2sub(size(theta), loca);
    
    prot(iProt,1) = R;
    prot(iProt,2) = C;
    
    %Update the current density in each bin
    j_bins(n_idx) = j_bins(n_idx) + max_val;
    
    count = zeros(n,1);
    for iCount = 1:size(prot,1)
        a = prot(iCount,1);
        b = prot(iCount,2);
        
        count(a) = count(a) + 1;
        count(b) = count(b) + 1;
        
        M(prot(iCount,1),prot(iCount,2)) = 0;
        M(prot(iCount,2),prot(iCount,1)) = 0;
    end
    
    A = find(count >= 2);
    
    M(A,:) = 0;
    M(:,A) = 0;
    
end

%%
%Find total current density over whole protocol - the larger the better
total_j = sum(j_bins);
%Find the total number of bins that don't have any current denisty going
%through them. Want this as small as possible
zero_bins = length(find(j_bins == 0));

%%
%Adjust the protocol so the larger value comes first, just because of the
%way w_mag, w_x, w_y and w_z is upper triangular matrix
prt=zeros(n-1,2);
for i= 1:n-1
    if prot(i,1) > prot(i,2)
        prt(i,1) = prot(i,2);
        prt(i,2) = prot(i,1);
    else
        prt(i,1) = prot(i,1);
        prt(i,2) = prot(i,2);
    end
end

%%
<<<<<<< HEAD:src/Matlab/Optimal Protocol Direction/angle_with_order.m
%Plot the vectors
figure
for i = 1:n-1
    quiver3(0,abs(w_x(prt(i,1),prt(i,2))), abs(w_y(prt(i,1),prt(i,2))), abs(w_z(prt(i,1),prt(i,2))));
    hold on;
end

%%
%Plot normalised vectors within unit sphere to make sure they are covering
%1/8th of area of sphere
figure
for i = 1:n-1
    quiver3(0,0,0,abs(w_x(prt(i,1),prt(i,2))/rho(prt(i,1),prt(i,2))), abs(w_y(prt(i,1),prt(i,2))/rho(prt(i,1),prt(i,2))), abs(w_z(prt(i,1),prt(i,2))/rho(prt(i,1),prt(i,2))));
    hold on;
=======
figure
for i = 3%:n-1
    if w_x(prt(i,1),prt(i,2)) < 0
quiver3(0,0,0,-1*(w_x(prt(i,1),prt(i,2))/rho(prt(i,1),prt(i,2))),-1*(w_y(prt(i,1),prt(i,2))/rho(prt(i,1),prt(i,2))),-1*(w_z(prt(i,1),prt(i,2))/rho(prt(i,1),prt(i,2))));
    else
        quiver3(0,0,0,(w_x(prt(i,1),prt(i,2))/rho(prt(i,1),prt(i,2))),(w_y(prt(i,1),prt(i,2))/rho(prt(i,1),prt(i,2))),(w_z(prt(i,1),prt(i,2))/rho(prt(i,1),prt(i,2))));
    end
hold on;
>>>>>>> master:src/Matlab/Optimal Protocol Direction/direction_hemisphere.m
end

[x, y, z] = sphere(128);
h = surfl(x, y, z);
set(h, 'FaceAlpha', 0.3)
shading interp

daspect([1 1 1]);
axis vis3d
set(gca,'XDir','rev','YDir','rev');

%%
%Plots a heat map on eighth sphere
%N.B have changed contents of sphere3D so the colour scale is correct
if ~isempty(rm)
    rm =sort(rm);
    m = fliplr([1:length(rm)]) - 1;
    j_temp = j_bins;
    j_idx = zeros(n_angles*n_angles,1);
    for iRm = 1:length(rm)
        j_idx(1:rm(iRm)-1) = j_temp(1:rm(iRm) -1);
        j_idx(rm(iRm)) = 0;
        j_idx(rm(iRm)+1:end-m(iRm)) = j_temp(rm(iRm):end);
        j_temp = j_idx(1:end-m(iRm));
    end
else
    j_idx = j_bins;
end

j_idx = reshape(j_idx, n_angles, n_angles);
<<<<<<< HEAD:src/Matlab/Optimal Protocol Direction/angle_with_order.m

=======
>>>>>>> master:src/Matlab/Optimal Protocol Direction/direction_hemisphere.m

%%
coord = zeros(n_angles+1,n_angles+1);
coord(1:end-1,1:end-1) = j_idx;
figure;
[xout,yout,zout,cout]=sphere3d(coord,0,pi,0,pi,1,1, 'surf',.0000000000001); %original .000001

cout = cout(2:end,1:end-1);

%%
figure;
for ix = 1:n_angles
    for iy = 1:n_angles
        c = [j_idx(iy,ix),j_idx(iy,ix);j_idx(iy,ix),j_idx(iy,ix)];
        surf(xout(ix:ix+1, iy:iy+1), yout(ix:ix+1, iy:iy+1), zout(ix:ix+1, iy:iy+1), c);
        hold on;
    end
end

daspect([1 1 1]);
axis vis3d
set(gca,'XDir','rev','YDir','rev');

