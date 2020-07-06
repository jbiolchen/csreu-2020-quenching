function [zeta, phigrid, fftphigrid, cgrid] = fsolve_phi(N, ell, zeta, eps, c_x, k_y, options)
% previously named solve_int2d
% calculate phi and FT(phi) using Newton's method, where desired solution u = phi - zeta

cgrid = c_x*(2/3).^(0:19);
k_xgrid = zeros(length(cgrid),1);
phigrid = zeros(N,length(cgrid));
fftphigrid = zeros(N,length(cgrid));

k_xinit = 1; % initial k_x; notice that k_x = -c_y*k_y/c_x with c_y = -c_x/k_y
phiinit = zeros(N,1); % initial function
phiextinit = [phiinit;k_xinit]; % initial data

% Newton's method
for j = 1:length(cgrid)
  f = @(phiext) int_eq2d(phiext,cgrid(j),k_y,ell,eps,N,zeta);
  phiextnew = fsolve(f,phiextinit, options);
  phiextinit = phiextnew;
  k_xgrid(j) = phiextinit(end);
  phigrid(:,j) = phiextinit(1:end-1)';
  fftphigrid(:,j) = abs(fft(phiextinit(1:end-1)))';
end
end