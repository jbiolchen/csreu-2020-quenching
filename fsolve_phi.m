function [zeta, rtn1, rtn2, cgrid] = fsolve_phi(N, L, ell, zeta, eps, c_x, k_y, options, regrid_error)
% regrid_error determines the behavior/output of fsolve_phi:
% regrid_error < 0 turns regridding off; phigrid, fftphigrid returned
% regrid_error >= 0 is the threshold for regridding (doubling N); k_xgrid returned

cgrid = [c_x*(6/7).^(0:39) 0];
k_xgrid = zeros(1, length(cgrid)); % row vector so other functions can define [k_xgrid; cgrid]
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
    if regrid_error < 0
        phigrid(:,j) = phiextinit(1:end-1)';
        fftphigrid(:,j) = abs(fft(phiextinit(1:end-1)))';
    else
        err=abs(fft(phiextinit(1:end-1)));
        err=max(err(3*N/8:5*N/8));
        if err > regrid_error
            display(['now refining grid to new N=' num2str(N*2)])
            phih=fft(phiextinit(1:end-1));
            phih = [phih(1:N/2);zeros(N,1);phih(1+N/2:N)];
            phiextinit = [2*ifft(phih,'symmetric');phiextinit(N+1)]; % REGRID
            N=2*N;
            ell = [[0:N/2] [-N/2+1: -1]]'; % variable in Fourier space
            zeta = (0:L/N:L-L/N)';
            f = @(phiext) int_eq2d(phiext,cgrid(j),k_y,ell,eps,N,zeta);
            phiextnew = fsolve(f,phiextinit, options);
            phiextinit = phiextnew;
            k_xgrid(j) = phiextinit(end);
        end   
    end
end

if regrid_error < 0
    rtn1 = phigrid; % desired solution u = phi - zeta
    rtn2 = fftphigrid;
else
    rtn1 = k_xgrid;
    rtn2 = NaN;
end
end