function plot_u(zeta, phigrid, cgrid)
% pass in zeta, phigrid, cgrid from fsolve_phi
% plot the desired solution u = phi - zeta

colors = distinguishable_colors(length(cgrid));

hold on
for m = 1:length(cgrid)
  plot(zeta,min(phigrid(:,m)-zeta,0), 'Color', colors(m, :))
end
hold off
title(legend(string(cgrid)), 'c_x')
legend("Location", "eastoutside")
xlabel('zeta')
ylabel('v(0,zeta)')
end
