function [type, ind] = seq_type(seq)
% given an array, return whether it's monotone increasing, monotone decreasing, or neither
% using 1, 1i, 0 for type lets you count the number of each type of sequence in a collection of arrays using sum, real, and imag
% if non-monotone, return the smallest index corresponding to the global extremum of interest

diffs = diff(seq);
ind = -1; % default value for monotone case (i.e. no such index)
if min(diffs) >= 0 % monotone increasing (including constant case)
    type = 1;
elseif max(diffs) <= 0 % monotone decreasing
    type = 1i;
else % non-monotone
    type = 0;
    for i = 1:length(diffs)
        if diffs(i) < 0 % concave up
            [~, ind] = min(seq);
            break
        elseif diffs(i) > 0 % concave down
            [~, ind] = max(seq);
            break
        else
            continue
        end
    end
end
end
