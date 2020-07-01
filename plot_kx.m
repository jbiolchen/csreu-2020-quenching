function plot_kx (k_ys, kx_array, kx_param)
% k_ys, kx_array are the output of fsolve_kx
% kx_param specifies what k_x is a function of, "c_x" or "k_y" (not yet needed)

colors = distinguishable_colors(length(k_ys));
hold on
for i = 1:length(k_ys)
    plot(kx_array(2*i,:),kx_array(2*i-1,:), 'Color', colors(i, :))
end
hold off
title(legend(string(k_ys)), 'k_y')
legend("Location", "eastoutside")
xlabel('c_x')
ylabel('k_x')
title('k_x vs. c_x for fixed k_y')

end