function [prot,prt] = max_EP_prot(nm, el_m, bad_chan, max_EP, map)

%nm - electrode choice as a string e.g 'LA' or 'RB'
%el_m - electrode choice as variable e.g LA or RA that are in 2_switch_map
%bad_chan - bad electrodes
%max_EP - electrode with maximum EP
%map - variable in 2_switch_map


load('pos.mat', 'pos');

%Recognize whether left or right electrode is being used and get
%corresponding positions
if nm == num2str('LA') | nm == num2str('LB')
    pos = pos(1:57,:);
    map_n = map(1:57,:);
    el = el_m;
else 
    pos = pos(58:end,:);
    map_n = map(58:end-1,:);
    el = el_m - 57;
end

%Just find coordinates of part of array we are interested in
h = pos(el,:);

%For plotting at the end so we can visualise whole array
n_el = setdiff([1:57],el)';
h_n = pos(n_el,:);

map_n = map_n(el);

%bad = bad_chan;
if bad_chan
bad_chan = bad_chan';
for i = 1:size(bad_chan,1)
    bad(i) = find(bad_chan(i) == map_n);
end
bad = bad'; 
else 
    bad = [];
end
    

m_EP = find(max_EP == map_n);

%For first iteration the n_id is the electrode with max EP
n_id = m_EP;

%Electrodes to ignore are bad electrodes
ig_el = bad;


j=1;
prt = [];
A = [];

%Calculate the distance from each electrode to the one with max EP
%Use to find go to each electrode
dist_mEP = (h-repmat(h(m_EP,:),length(h),1));
dist_mEP(ig_el,:) = 0;

for i = 1:size(dist_mEP,1)
    dist_abs_mEP(i) = norm(dist_mEP(i,:),2);
end
dist_abs_mEP = dist_abs_mEP';



while sum(dist_abs_mEP) ~= 0

dist = (h-repmat(h(n_id,:),length(h),1));
dist(ig_el,:) = 0;

for i = 1:size(dist,1)
    dist_abs(i) = norm(dist(i,:),2);
end
dist_abs = dist_abs';

%Sort electrodes in terms of largest distance from electrode of interest
[max_d, max_i] = sort(dist_abs,'descend');

%Take the third furthest electrode
prt(j,1) = n_id;
prt(j,2) = max_i(3);

j=j+1;


%Find all electrodes that are not adjacent to the electrode of interest
%Adjacent electrodes will be seperated by 1.3mm
%N.B. Also automatically ignores all ig_el as these have dist_abs = 0
id = find(dist_abs >1.31);

%Sort electrodes in terms of smallest distance from electrode of interest
[min_d, min_i] = sort(dist_abs(id));

if size(min_i,1) < 2
    break
end

%Take the second nearest electrode
min_i = id(min_i(2));
prt(j,1) = n_id;
prt(j,2) = min_i;

j=j+1;

%Calculate how many times an electrode has been used in the protocol so far
elec = zeros(length(el),1);
for i = 1:size(prt,1)
    a = prt(i,1);
    b = prt(i,2);
    
    elec(a) = elec(a) + 1;
    elec(b) = elec(b) + 1;
end

%If an electrode has been used four times, no longer want to inject through
%this electrode so add to list of electrodes to ignore for next iteration
A = find(elec >= 2);
ig_el = [bad;A];

%Set to zero so they are not considered when finding the next electrode to
%start from
dist_abs_mEP(ig_el) = 0;
%Set the electrode you were just considering to zero too so it doesnt go
%back and forth
dist_abs_mEP(n_id) = 0;

%Need to find the electrode that is next closest to the m_EP
%Find the index with the minimum non zero distance
id2 = find(dist_abs_mEP>0);
[min_id,n_id] = min(dist_abs_mEP(id2));
n_id = id2(n_id);

%As transpose mucks things up so best to clear each time
clear dist_abs

end

prot = map(el_m(prt));

figure
for i = 1:size(prt,1)
subplot(6,6,i);
scatter(h(:,1),h(:,2), 'b', 'MarkerFaceColor', 'b');
hold on;
scatter(h_n(:,1),h_n(:,2), 'b');
hold on;
scatter(h(prt(i,:),1), h(prt(i,:),2), 'r', 'MarkerFaceColor', 'r');
hold on;
scatter(h(m_EP,1),h(m_EP,2), 'g', 'MarkerFaceColor', 'g');
end


    
