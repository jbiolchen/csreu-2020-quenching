function [const_id, constlist, var_id, kx_array] = fsolve_kx(const_id, constlist, var_id, var, N, numiter, ell, zeta, eps, k_xinit, phiinit, ds, dvar, options)
% const is k_y or c_x; var is the other one
% for each const in constlist, calculate k_x as a function of var using secant continuation

const_ids = ["k_y", "c_x"]; % possible values of const_id

kx_array = zeros(2*length(constlist), numiter);
for i = 1:length(constlist)
    const = constlist(i);
    
    % "baby continuation" to calculate first secant for secant continuation
    % calculate first point
    f_map = containers.Map(const_ids, {@(phiext) int_eq2d(phiext,var,const,ell,eps,N,zeta), @(phiext) int_eq2d(phiext,const,var,ell,eps,N,zeta)});
    f = f_map(const_id);
    phiextinit = [phiinit;k_xinit]; % initial data
    phiextsol = fsolve(f, phiextinit, options);
    uold = [phiextsol; var];
    var=var+dvar;
    f_map = containers.Map(const_ids, {@(phiext) int_eq2d(phiext,var,const,ell,eps,N,zeta), @(phiext) int_eq2d(phiext,const,var,ell,eps,N,zeta)});
    f = f_map(const_id);
    
    % calculate second point
    phiextsol = fsolve(f, phiextsol, options);
    u0 = [phiextsol; var];

    % calculate first secant
    sec = u0 - uold;
    sec = sec/norm(sec);

    % secant continuation
    for j = 1:numiter
        uinit = u0 + sec*ds;
        fext = @(u) int_eq2d_ext(u,const_id,const,ell,eps,N,zeta,sec,uinit);   
        usol = fsolve(fext, uinit, options);
        kx_array(2*i-1:2*i, j) = usol(end-1:end);
        sec = usol - u0;
        sec = sec/norm(sec);
        u0=usol;
    end
end
end