function f = int_eq2d_ext(phisec,k_y,ell,eps,N,zeta,sec,uinit)
% used in fsolve_kx
% applies boundary integral operator to phisec and concatenates dot product

c_x = phisec(end);
phiext = phisec(1:end-1);
f0 = int_eq2d(phiext,c_x,k_y,ell,eps,N,zeta);
f1 = dot(phisec-uinit,sec);
f = [f0;f1];
end
