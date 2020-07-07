function plot_extrloc(constlist, kx_array, extrema)
% constlist, kx_array are outputs of fsolve_kx
% extrema is an output of seqtype_kx
% for non-monotone plots of k_x vs. c_x, plots (k_y, c_x) for each global extremum of interest

color = distinguishable_colors(1);
[num, ~] = size(extrema);
xlist = zeros(num, 1);
ylist = zeros(num, 1);
for i = 1:num
    xlist(i, 1) = constlist(extrema(i, 1));
    ylist(i, 1) = kx_array(2*extrema(i, 1), extrema(i, 2));
end
plot(xlist, ylist, 'Color', color)
xlabel('k_y')
ylabel('c_x')
title("location of global minima of non-monotone k_x(c_x)")
end