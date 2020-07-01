function plot_kx (const_id, constlist, var_id, kx_array)
% [const_id, constlist, var_id, kx_array] is the output of fsolve_kx

colors = distinguishable_colors(length(constlist));
hold on
for i = 1:length(constlist)
    plot(kx_array(2*i,:),kx_array(2*i-1,:), 'Color', colors(i, :))
end
hold off
title(legend(string(constlist)), const_id)
legend("Location", "eastoutside")
xlabel(var_id)
ylabel('k_x')
title("k_x vs. " + var_id + " for fixed " + const_id)

end