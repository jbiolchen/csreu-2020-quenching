function type = seq_type(seq)
% given an array, return whether it's monotone increasing, monotone decreasing, or neither
% returning 1, 1i, 0 lets you count the number of each type of sequence in a collection of arrays using sum, real, and imag

diffs = diff(seq);
if min(diffs) >= 0 % monotone increasing
    type = 1;
elseif max(diffs) <= 0 % monotone decreasing
    type = 1i;
else
    type = 0; % not monotone
end
end