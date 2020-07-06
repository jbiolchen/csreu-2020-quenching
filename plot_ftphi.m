function plot_ftphi(zeta, fftphigrid, cgrid, N)
% pass in zeta, fftphigrid, cgrid from fsolve_phi
% plot the Fourier transform of phi

colors = distinguishable_colors(length(cgrid));

hold on
for m = 1:length(cgrid)
  %plot(zeta,fftphigrid(:,m), 'Color', colors(m, :))
  plot(zeta(N/4:3*N/4),fftphigrid(N/4:3*N/4,m), 'Color', colors(m, :))
end
hold off
title(legend(string(cgrid)), 'c')
legend("Location", "eastoutside")
xlabel('zeta')
ylabel('fft(phi)')
% xlim([pi/2 3*pi/2])
end