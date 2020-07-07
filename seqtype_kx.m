function [const_inds, kx_types, extrema] = seqtype_kx(given_type, kx_array, constlist)
% kx_array, constlist are outputs of fsolve_kx
% given_type is one of the possible first outputs of seq_type
% const_inds contains the indices of constlist corresponding to kx_array curves of type given_type
% kx_types contains the curve types for all curves in kx_array
% extrema contains the (const, k_x) indices of the global extrema of interest

numconsts = length(constlist);
kx_types = zeros(numconsts, 1);
const_inds = [];
extrema = []; % will remain empty if given_type == 1 or 1i
for i = 1:numconsts
    [type, extr_ind] = seq_type(kx_array(2*i-1,:));
    kx_types(i, 1) = type;
    if type == given_type
        const_inds = [const_inds; i];
        if given_type == 0
            extrema = [extrema; i extr_ind];
        end
    end     
end
end