function [const_inds, kx_types] = seqtype_kx(given_type, kx_array, constlist)
% kx_array, constlist are outputs of fsolve_kx
% given_type is one of the possible outputs of seq_type
% return the indices of constlist corresponding to kx_array curves of type given_type, and a list of all curve types

numconsts = length(constlist);
kx_types = zeros(numconsts, 1);
const_inds = [];
for i = 1:numconsts
    type = seq_type(kx_array(2*i-1,:));
    kx_types(i, 1) = type;
    if type == given_type
        const_inds = [const_inds; i];
    end     
end
end