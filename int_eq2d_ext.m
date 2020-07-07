function f = int_eq2d_ext(phisec,const_id,const,ell,eps,N,zeta,sec,uinit)
% used in fsolve_kx
% applies boundary integral operator to phisec and concatenates dot product

const_ids = ["k_y", "c_x"]; % possible values of const_id

var = phisec(end);
phiext = phisec(1:end-1);
f0_map = containers.Map(const_ids, {int_eq2d(phiext,var,const,ell,eps,N,zeta), int_eq2d(phiext,const,var,ell,eps,N,zeta)});
f0 = f0_map(const_id);
f1 = dot(phisec-uinit,sec);
f = [f0;f1];
end
