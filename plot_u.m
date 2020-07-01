function plot_u(zeta, phigrid, cgrid)
% pass in zeta, phigrid, cgrid from fsolve_phi
% plot the desired solution u = phi - zeta

colors = distinguishable_colors(length(cgrid));

hold on
for m = 1:length(cgrid)
  plot(zeta,phigrid(:,m)-zeta, 'Color', colors(m, :))
end
hold off
title(legend(string(cgrid)), 'c')
legend("Location", "eastoutside")
xlabel('zeta')
ylabel('u')
end