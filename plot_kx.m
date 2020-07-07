function plot_kx (const_id, constlist, var_id, kx_array, const_inds, const_inds2, points)
% [const_id, constlist, var_id, kx_array] is the output of fsolve_kx
% const_inds, const_inds2 specify kx_array curves to be plotted with different line styles, and should have no common elements
% each row of points is a (const, k_x) index pair to be plotted with a marker

colors = distinguishable_colors(length(constlist));
hold on
for i = 1:length(constlist)
    % plot curve i
    if ismember(i, const_inds)
        plot(kx_array(2*i,:),kx_array(2*i-1,:), '--', 'Color', colors(i, :)) % dashed line
    elseif ismember(i, const_inds2)
        plot(kx_array(2*i,:),kx_array(2*i-1,:), '-.', 'Color', colors(i, :)) % dash-dotted line
    else
        plot(kx_array(2*i,:),kx_array(2*i-1,:), 'Color', colors(i, :)) % solid line
    end
    
    % check for a point on curve i that should be plotted with a marker
    if ismember(i, points(:, 1))
        col = points(find(points(:,1) == i), 2);
        plot(kx_array(2*i, col), kx_array(2*i-1, col), 'o', 'Color', colors(i, :), 'HandleVisibility', 'off')
    end  
end
hold off
title(legend(string(constlist)), const_id)
legend("Location", "eastoutside")
xlabel(var_id)
ylabel('k_x')
title("k_x vs. " + var_id + " for fixed " + const_id)
end