% main function where constants are defined and other functions are called

close all;
clear all;
 
% define optimization parameters
N = 2^9; %  number of Fourier modes
numiter = 100; % number of secant continuation iterations
ds = 5/sqrt(N); % secant continuation step size
dvar = -1e-4; % "baby continuation" step size/direction
options = optimset('Jacobian','off','Display','iter','TolFun',1e-6,'TolX',1e-6,'MaxIter',10,'Algorithm','trust-region-reflective');
regrid_error = 1e-4; % threshold to double N in fsolve_phi, fsolve_kx
adaptive_ds = true;
jac_error = 1e-6;

% define constants
eps = 0.3; % theta parameter
c_x = 0.1;
c_xs = [0, 0.00001, 0.0001, 0.001, 0.01]; % fixed values for fsolve_kx only
k_y = 1;
k_ys = 0.1:0.1:0.6; % fixed values for fsolve_kx only
L = 2*pi; % length of domain
k_xinit = 1; % initial k_x; notice that k_x = -c_y*k_y/c_x with c_y = -c_x/k_y
phiinit = zeros(N,1); % initial function
ell = [0:N/2 -N/2+1:-1]'; % variable in Fourier space
zeta = (0:L/N:L-L/N)';


% ========== SAMPLE CODE FOR USING ALL OF THE OTHER FUNCTIONS ==========
% % plot u = phi - zeta at x = 0 (with regridding off)
% [fsp1, fsp2, fsp3, fsp4] = fsolve_phi(N, L, ell, zeta, eps, c_x, k_y, options, -1);
% figure('Name', 'plot_u')
% plot_u(fsp1, fsp2, fsp4);
% 
% % plot FT(phi) vs. zeta
% figure('Name', 'plot_ftphi')
% plot_ftphi(fsp1, fsp3, fsp4);
% 
% % plot k_x vs. c_x for each k_y in k_ys
% [fsk1, fsk2, fsk3, fsk4] = fsolve_kx("k_y", k_ys, "c_x", c_x, N, L, numiter, ell, zeta, eps, k_xinit, phiinit, ds, dvar, options, regrid_error, adaptive_ds, jac_error);
% figure('Name', 'plot_kx in c_x')
% plot_kx(fsk1, fsk2, fsk3, fsk4, [], [], [])
%
% % as above, with monotone plots and global minima of non-monotone plots marked
% monoinc = seqtype_kx(1, fsk4, fsk2);
% monodec = seqtype_kx(1i, fsk4, fsk2);
% [~, ~, pts] = seqtype_kx(0, fsk4, fsk2);
% figure('Name', 'plot_kx in c_x, monotone marked')
% plot_kx(fsk1, fsk2, fsk3, fsk4, monoinc, monodec, pts)
% 
% % for non-monotone plots of k_x vs. c_x, plot the (k_y, c_x) corresponding to global minima
% figure('Name', 'plot_extrloc')
% plot_extrloc(fsk2, fsk4, pts)
% 
% % plot k_x vs. k_y for each c_x in c_xs
% [fsk5, fsk6, fsk7, fsk8] = fsolve_kx("c_x", c_xs, "k_y", k_y, N, L, numiter, ell, zeta, eps, k_xinit, phiinit, ds, dvar, options, regrid_error, adaptive_ds, jac_error);
% figure('Name', 'plot_kx in k_y')
% plot_kx(fsk5, fsk6, fsk7, fsk8, [], [], [])
% 
% % given k_y, plot k_x vs. c_x with regridding on
% % this functionality is better implemented in fsolve_kx
% [zeta, k_xgrid, ~, cgrid] = fsolve_phi(N, L, ell, zeta, eps, c_x, k_y, options, regrid_error);
% figure('Name', "fsolve_phi: k_x vs. c_x for k_y = " + string(k_y))
% plot_kx("k_y", [k_y], "c_x", [k_xgrid; cgrid], [], [], [])
% 
% % plot H^{1/2} (Sobolev) norm vs. k_y for c_x = 0
% [~, ~, ~, ~, norm_array] = fsolve_kx("c_x", [0], "k_y", k_y, N, L, numiter, ell, zeta, eps, k_xinit, phiinit, ds, dvar, options, regrid_error, adaptive_ds, jac_error);
% figure
% plot(norm_array(:, 1), norm_array(:, 2), '-o', 'Color', [1 0 0])
% xlabel('k_y')
% ylabel('norm')
% title('H^1^/^2 norm vs. k_y for c_x = 0')
% ======================================================================

% ===================== SAMPLE CODE FOR DEBUGGING ======================
% % check for errors in jac_inteq2d Jacobian formulas
% % (correct Jacobian produces err quadratic in du)
% const=0;
% sec=rand(N+2,1)*.01;
% u=cos(ones(N+2,1));
% du=1e-5*rand(N+2,1);
% % du(end-1:end)=0; % for isolating Jacobian errors
% F0 = @(phisec) int_eq2d_ext(phisec,"k_y",const,ell,eps,N,zeta,sec,u);
% J = @(phisec, dphisec) jac_inteq2dext(dphisec,phisec,"k_y",const,ell,eps,N,zeta,sec);
% err = F0(u+du)-F0(u)-J(u,du);
% plot(err)
% ======================================================================


% various Sobolev norm plots
[~, ~, ~, ~, norm_array] = fsolve_kx("c_x", [0], "k_y", k_y, N, L, numiter, ell, zeta, eps, k_xinit, phiinit, ds, dvar, options, regrid_error, adaptive_ds, jac_error);
figure
plot(norm_array(:, 1), norm_array(:, 2), '-o', 'Color', [1 0 0])
xlabel('k_y')
ylabel('norm')
title('H^1^/^2 norm vs. k_y for c_x = 0')

figure
plot(log(norm_array(:, 1)), log(norm_array(:, 2)), '-o', 'Color', [1 0 0])
xlabel('ln(k_y)')
ylabel('ln(norm)')
title('ln(H^1^/^2 norm) vs. ln(k_y) for c_x = 0')

figure
plot(norm_array(:, 1), norm_array(:, 2)./norm_array(:, 1), '-o', 'Color', [1 0 0])
xlabel('k_y')
ylabel('norm')
title('1/k_y * H^1^/^2 norm vs. k_y for c_x = 0')

figure
plot(log(norm_array(:, 1)), log(norm_array(:, 2)./norm_array(:, 1)), '-o', 'Color', [1 0 0])
xlabel('ln(k_y)')
ylabel('ln(norm)')
title('ln(1/k_y * H^1^/^2 norm) vs. ln(k_y) for c_x = 0')

% plot k_x vs. k_y for each c_x in c_xs
% dotted line = monotone increasing, dot-dashed line = monotone decreasing
[const_id, constlist, var_id, kx_array] = fsolve_kx("c_x", c_xs, "k_y", k_y, N, L, numiter, ell, zeta, eps, k_xinit, phiinit, ds, dvar, options, regrid_error, adaptive_ds, jac_error);
monoinc = seqtype_kx(1, kx_array, constlist);
monodec = seqtype_kx(1i, kx_array, constlist);
[~, ~, extr_pts] = seqtype_kx(0, kx_array, constlist);
figure('Name', 'plot_kx in c_x')
plot_kx(const_id, constlist, var_id, kx_array, monoinc, monodec, extr_pts)
