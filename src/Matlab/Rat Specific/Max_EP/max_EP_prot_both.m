function [prt, prot] = max_EP_prot_both(nm, bad_chan_A, bad_chan_B, max_EP, max_EP_elec)

%nm - electrode side you want to use as a string e.g 'L' or 'R'
%bad_chan_A - bad electrodes on A side from 1-32
%bad_chan_B - bad electrodes on B side from 1-32
%max_EP - index of electrode with max EP
%max_EP_elec - side max EP is on, e.g LA or LB

load('pos.mat', 'pos');
load('2_switch_map.mat');


%Recognize whether left or right electrode is being used and get
%corresponding positions
if nm == num2str('L') 
    pos = pos(1:57,:);
    A = LA;
    B = LB;
else 
    pos = pos(58:end,:);
    A = RA;
    B = RB;
end

map_n = map;
h = pos;

%For plotting at the end so we can visualise whole array
%n_el = setdiff([1:57],el)';
%h_n = pos(n_el,:);

%map_n = map_n(el);

%bad = bad_chan;
if bad_chan_A
bad_chan_A = bad_chan_A';
    for i = 1:size(bad_chan_A,1)
        bad = find(bad_chan_A(i) == map_n);
        bad_A(i) = bad(find(ismember(bad,A)));
    end
bad_A = bad_A'; 
else
    bad_A = [];
end

if bad_chan_B
bad_chan_B = bad_chan_B';
    for i = 1:size(bad_chan_B,1)
        bad = find(bad_chan_B(i) == map_n);
        bad_B(i) = bad(find(ismember(bad,B)));
    end
bad_B = bad_B'; 
else
    bad_B = [];
end

bad = [bad_A;bad_B];


m_EP = find(max_EP == map_n);
id_EP = find(ismember(m_EP,max_EP_elec));
m_EP = m_EP(id_EP);

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
    sign_x(i) = sign(dist_mEP(i,1));
    sign_y(i) = sign(dist_mEP(i,2));
end
dist_abs_mEP = dist_abs_mEP';
sign_x = sign_x';
sign_y = sign_y';



while sum(dist_abs_mEP) ~= 0

dist = (h-repmat(h(n_id,:),length(h),1));
dist(ig_el,:) = 0;

for i = 1:size(dist,1)
    dist_abs(i) = norm(dist(i,:),2);
end
dist_abs = dist_abs';

sign_n_id = sign_x(n_id);


%Sort electrodes in terms of largest distance from electrode of interest
[max_d, max_i] = sort(dist_abs,'descend');
sign_max_i = sign_x(max_i);
op_sign= find(sign_max_i == -1*sign_n_id);


if op_sign
prt(j,1) = n_id;
prt(j,2) = max_i(op_sign(2));
else
%Take the third furthest electrode
prt(j,1) = n_id;
prt(j,2) = max_i(4);
 end

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

sign_min_i = sign_x(min_i);
op_sign = find(sign_min_i == -1*sign_n_id);

if op_sign 
min_i = id(min_i(op_sign(1)));
prt(j,1) = n_id;
prt(j,2) = min_i;
else

%Take the second nearest electrode
min_i = id(min_i(2));
prt(j,1) = n_id;
prt(j,2) = min_i;
 end

j=j+1;

%Calculate how many times an electrode has been used in the protocol so far
elec = zeros(length(pos),1);
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

prot = map_n(prt);


good = setdiff([1:57],bad);
h_g = pos(good,:);
h_n = pos(bad,:);

figure
for i = 1:size(prt,1)
%for i = 1:36
subplot(6,12,i);
scatter(h_g(:,1),h_g(:,2), 'b', 'MarkerFaceColor', 'b');
hold on;
scatter(h_n(:,1),h_n(:,2), 'b');
hold on;
scatter(h(prt(i,:),1), h(prt(i,:),2), 'r', 'MarkerFaceColor', 'r');
hold on;
scatter(h(m_EP,1),h(m_EP,2), 'g', 'MarkerFaceColor', 'g');
end


    
