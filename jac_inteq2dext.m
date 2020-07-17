function Df = jac_inteq2dext(dphisec,phisec,const,ell,eps,N,zeta,sec)
% for adding Jacobian to fsolve

var = phisec(end);%cx
phiext = phisec(1:end-1);

dvar = dphisec(end);% dcx
dphiext = dphisec(1:end-1);

Df0 = jac_inteq2d(dphiext,phiext,var,dvar,const,ell,eps,N,zeta);
Df1=dot(dphisec,sec);
Df = [Df0;Df1];
end