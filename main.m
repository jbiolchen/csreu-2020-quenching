% main function where constants are defined and other functions are called

close all;
clear all;

% define optimization parameters
N = 2^9; %  number of Fourier modes
numiter = 30; % number of secant continuation iterations
ds = 0.8; % secant continuation step size
dvar = 1e-4; % "baby continuation" step size/direction
options = optimset('Jacobian','off','Display','iter','TolFun',1e-6,'TolX',1e-6,'MaxIter',10,'Algorithm','trust-region-reflective');

% define constants
eps = 0.3; % theta parameter
c_x = 0.1;
c_xs = 0.1:0.1:1.1; % fixed values for fsolve_kx only
k_y = 0.1;
k_ys = 0.1:0.1:1.1; % fixed values for fsolve_kx only
L = 2*pi; % length of domain
k_xinit = 1; % initial k_x; notice that k_x = -c_y*k_y/c_x with c_y = -c_x/k_y
phiinit = zeros(N,1); % initial function
ell = [[0:N/2] [-N/2+1: -1]]'; % variable in Fourier space
zeta = (0:L/N:L-L/N)';


% ========== sample code for using all of the other functions ==========
% % plot u = phi - zeta at x = 0
% [fsp1, fsp2, fsp3, fsp4] = fsolve_phi(N, ell, zeta, eps, c_x, k_y, options);
% figure('Name', 'plot_u')
% plot_u(fsp1, fsp2, fsp4);
% 
% % plot FT(phi) vs. zeta
% figure('Name', 'plot_ftphi')
% plot_ftphi(fsp1, fsp3, fsp4);
% 
% % plot k_x vs. c_x for each k_y in k_ys
% [fsk1, fsk2, fsk3, fsk4] = fsolve_kx("k_y", k_ys, "c_x", c_x, N, numiter, ell, zeta, eps, k_xinit, phiinit, ds, dvar, options);
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
% [fsk5, fsk6, fsk7, fsk8] = fsolve_kx("c_x", c_xs, "k_y", k_y, N, numiter, ell, zeta, eps, k_xinit, phiinit, ds, dvar, options);
% figure('Name', 'plot_kx in k_y')
% plot_kx(fsk5, fsk6, fsk7, fsk8, [], [], [])
% ======================================================================


% plot k_x vs. c_x for each k_y in k_ys
% dotted line = monotone increasing, dot-dashed line = monotone decreasing
[const_id, constlist, var_id, kx_array] = fsolve_kx("k_y", k_ys, "c_x", c_x, N, numiter, ell, zeta, eps, k_xinit, phiinit, ds, dvar, options);
monoinc = seqtype_kx(1, kx_array, constlist);
monodec = seqtype_kx(1i, kx_array, constlist);
[~, ~, extr_pts] = seqtype_kx(0, kx_array, constlist);
figure('Name', 'plot_kx in c_x')
plot_kx(const_id, constlist, var_id, kx_array, monoinc, monodec, extr_pts)

% for non-monotone plots of k_x vs. c_x, plot the (k_y, c_x) corresponding to global minima
figure('Name', 'plot_extrloc')
plot_extrloc(constlist, kx_array, extr_pts)

