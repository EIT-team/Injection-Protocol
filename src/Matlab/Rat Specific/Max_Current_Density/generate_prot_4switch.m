function [Prot, Prot_ext, good_chan,prt] = generate_prot_4switch(bad_LA, bad_LB, bad_RA, bad_RB, w, n)

load('4_switch_map', 'fwd_map', 'inv_map');

bad_LB = bad_LB + 32;
bad_RA = bad_RA + 64;
bad_RB = bad_RB + 96;

bad = [bad_LA, bad_LB, bad_RA, bad_RB];
bad_chan = inv_map(bad);

elec = [1:114]';
good_chan = setdiff(elec, bad_chan);

w_elec = w(good_chan, good_chan);

[prot_max] = max_protocol(w_elec);
[prt] = ind_max_prot(prot_max, n);

prt = good_chan(prt);
Prot = fwd_map(prt);

%To get a few extra that give the maximum current density in region of
%interest
elec = [1:115]';
good_chan = setdiff(elec, bad_chan);
w_elec = w(good_chan, good_chan);
[prot_max] = max_protocol(w_elec);
Prot_ext = prot_max(1:20,:);
Prot_ext = good_chan(Prot_ext);
Prot_ext = fwd_map(Prot_ext);

load('pos.mat', 'pos');

for i = 1:size(prt,1)
    subplot(10,12,i)
    scatter(pos(:,1), pos(:,2), 'b', 'MarkerFaceColor', 'b');
    hold on;
    scatter(pos(prt(i,:),1), pos(prt(i,:),2), 'r', 'MarkerFaceColor', 'r');
end


% for i = 1:size(Prot,1)
% if Prot(i,1) > 64
% Prot(i,1) = Prot(i,1) - 32;
% end
% if Prot(i,2) > 64
% Prot(i,2) = Prot(i,2) - 32;
% end
% end