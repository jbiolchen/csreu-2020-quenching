% main function where constants are defined and other functions are called

close all;
clear all;

% define optimization parameters
N = 2^9; %  number of Fourier modes
numiter = 10; % number of secant continuation iterations
ds = 0.8; % secant continuation step size
dcx = 1e-4; % "baby continuation" step size/direction

% define constants
eps = 0.3; % theta parameter
c_x = 0.1;
k_y = 0.2;
L = 2*pi; % length of domain
k_xinit = 1; % initial k_x; notice that k_x = -c_y*k_y/c_x with c_y = -c_x/k_y
phiinit = zeros(N,1); % initial function
ell = [[0:N/2] [-N/2+1: -1]]'; % variable in Fourier space
zeta = (0:L/N:L-L/N)';


% sample code for using all of the other functions
% plot u = phi - zeta at x = 0
[fsp1, fsp2, fsp3, fsp4] = fsolve_phi(N, ell, zeta, eps, c_x, k_y);
figure('Name', 'plot_u')
plot_u(fsp1, fsp2, fsp4);

% plot FT(phi) vs. zeta
figure('Name', 'plot_ftphi')
plot_ftphi(fsp1, fsp3, fsp4);

% plot k_x vs. c_x for each k_y in k_ys
k_ys = 0.1:0.1:1;
[fsk1, fsk2] = fsolve_kx(N, numiter, ell, zeta, eps, c_x, k_ys, k_xinit, phiinit, ds, dcx);
figure('Name', 'plot_kx')
plot_kx(fsk1, fsk2, 'c_x')

