function [Prot, Prt, Prot_add, Prt_add, prt] = generate_prot_2switch(elec_L, elec_R, bad_chan_L, bad_chan_R, w, map, n)
%Gives out an independent protocol using each electrode n times
%
% input:    elec    choice of A or B set of electrodes on right and left
%                   electrode
%           bad_chan_L  index of bad channels on left electrode
%           bad_chan_R  index of bad channels in right electrode
%           w           matrix with all weights for different injection
%                       pairs for ROI
%           map         array that maps the electrode numbering from
%                       connectors to matlab simulations
%           n           the number of times you want to inject through each
%                       electrode( for independent must be 2)fwd_mapCopy


%% Reads in the bad channels in the electrode map format and converts them into the matlab map format
if elec_L 
    if elec_R
    elec = [elec_L;elec_R];
    else
        elec = elec_L;
    end
else
    elec = elec_R;
end


bad_chan = map_bad_chan(bad_chan_L, bad_chan_R, elec, map);

%% Finds the good channels and only uses these rows and columns in matrix w
good_chan = setdiff(elec,bad_chan);
w_elec = w(good_chan, good_chan);

%Lists the protocol lines in terms of maximum value of w
[prot_max] = max_protocol(w_elec);

%Finds the protcol that uses each electrode twice
[prt] = ind_max_prot(prot_max, n);
%[prt] = create_protocol(w_elec,56);

prt = good_chan(prt);
%prot_max = good_chan(prot_max);
Prot = map(prt);

for i = 1:size(prt,1)
    if prt(i,1) > 57 && prt(i,2) > 57
        Prt(i,1) = map(prt(i,1)) + 32;
        Prt(i,2) = map(prt(i,2)) + 32;
        
    elseif prt(i,1) <= 57 && prt(i,2) <= 57
        Prt(i,1) = map(prt(i,1));
        Prt(i,2) = map(prt(i,2));
        
    elseif prt(i,1) <=57 && prt(i,2) > 57
        Prt(i,1) = map(prt(i,1));
        Prt(i,2) = map(prt(i,2)) + 32;
        
    elseif prt(i,1) > 57 && prt(i,2) <= 57
        Prt(i,1) = map(prt(i,1)) + 32;
        Prt(i,2) = map(prt(i,2));
    end
end

% 
load('pos.mat', 'pos');
h = pos;
h_g = h(elec,:);
n_el = (setdiff([1:114], elec));
h_n = h(n_el,:);
h_b = h(bad_chan,:);
figure
for i = 1:size(prt,1)
%for i = 1:36
subplot(6,12,i);
scatter(h_g(:,1),h_g(:,2), 'b', 'MarkerFaceColor', 'b');
hold on;
scatter(h_n(:,1),h_n(:,2), 'b');
hold on;
scatter(h_b(:,1), h_b(:,2),'g');
hold on;
scatter(h(prt(i,:),1), h(prt(i,:),2), 'r', 'MarkerFaceColor', 'r');
%hold on;
%scatter(h(m_EP,1),h(m_EP,2), 'g', 'MarkerFaceColor', 'g');
end

elec = [elec; 115];
good_chan = setdiff(elec,bad_chan);
w_elec = w(good_chan, good_chan);
prot_max = max_protocol(w_elec);
prt_max = good_chan(prot_max);
Prot_add = map(prt_max);


for i = 1:size(prt_max,1)
    if prt_max(i,1) > 57 && prt_max(i,2) > 57
        Prt_add(i,1) = map(prt_max(i,1)) + 32;
        Prt_add(i,2) = map(prt_max(i,2)) + 32;
        
    elseif prt_max(i,1) <= 57 && prt_max(i,2) <= 57
        Prt_add(i,1) = map(prt_max(i,1));
        Prt_add(i,2) = map(prt_max(i,2));
        
    elseif prt_max(i,1) <=57 && prt_max(i,2) > 57
        Prt_add(i,1) = map(prt_max(i,1));
        Prt_add(i,2) = map(prt_max(i,2)) + 32;
        
    elseif prt_max(i,1) > 57 && prt_max(i,2) <= 57
        Prt_add(i,1) = map(prt_max(i,1)) + 32;
        Prt_add(i,2) = map(prt_max(i,2));
    
       
    end
end


