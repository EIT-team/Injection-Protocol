function [bad_chan] = map_bad_chan(bad_chan_L, bad_chan_R, elec, map)

D = [1:114]';
set_zero = setdiff(D, elec);
map_alt = map;
map_alt(set_zero) = 0;

L_elec = map_alt(1:57);
R_elec = map_alt(58:end);

bad_chan_L = bad_chan_L';
bad_chan_R = bad_chan_R';

if bad_chan_L 
    for i = 1:size(bad_chan_L)
        bc_L(i) = find(bad_chan_L(i) == L_elec);
    end
    
else bc_L = [];
end
    
if bad_chan_R
    for i = 1:size(bad_chan_R)
        bc_R(i) = find(bad_chan_R(i) == R_elec) + 57;
    end
  
else bc_R = [];
end

bad_chan = [bc_L,bc_R]';

