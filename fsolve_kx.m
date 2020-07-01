function [k_ys, kx_array] = fsolve_kx(N, numiter, ell, zeta, eps, c_x, k_ys, k_xinit, phiinit, ds, dcx)
% for each k_y in k_ys, calculate k_x as a function of c_x using secant continuation
% coming soon: given c_x, calculate k_x as a function of k_y using secant continuation

options = optimset('Jacobian','off','Display','iter','TolFun',1e-6,'TolX',1e-6,'MaxIter',10,'Algorithm','trust-region-reflective');

kx_array = zeros(2*length(k_ys),numiter);
for i = 1:length(k_ys)
    k_y = k_ys(i);

    % "baby continuation" to calculate first secant for secant continuation
    % calculate first point
    f = @(phiext) int_eq2d(phiext,c_x,k_y,ell,eps,N,zeta);
    phiextinit = [phiinit;k_xinit]; % initial data
    phiextsol = fsolve(f, phiextinit, options);
    uold = [phiextsol; c_x];
    c_x=c_x+dcx;
    f = @(phiext) int_eq2d(phiext,c_x,k_y,ell,eps,N,zeta);

    % calculate second point
    phiextsol = fsolve(f, phiextsol);
    u0 = [phiextsol; c_x];

    % calculate first secant
    sec = u0 - uold;
    sec = sec/norm(sec);

    % secant continuation
    for j = 1:numiter
        uinit = u0 + sec*ds;
        fext = @(u) int_eq2d_ext(u,k_y,ell,eps,N,zeta,sec,uinit);
        usol = fsolve(fext, uinit, options);
        kx_array(2*i-1:2*i, j) = usol(end-1:end);
        sec = usol - u0;
        sec = sec/norm(sec);
        u0=usol;
    end
end
end