function [const_id, constlist, var_id, kx_array] = fsolve_kx(const_id, constlist, var_id, var, N, L, numiter, ell, zeta, eps, k_xinit, phiinit, ds, dvar, options, regrid_error, adaptive_ds)
% const is k_y or c_x; var is the other one
% for each const in constlist, calculate k_x as a function of var using secant continuation

const_ids = ["k_y", "c_x"]; % possible values of const_id
phiextinit = [phiinit;k_xinit]; % initial data

kx_array = zeros(2*length(constlist), numiter);
for i = 1:length(constlist)
    const = constlist(i);
    
    % "baby continuation" to calculate first secant for secant continuation
    % calculate first point
    f_map = containers.Map(const_ids, {@(phiext) int_eq2d(phiext,var,const,ell,eps,N,zeta), @(phiext) int_eq2d(phiext,const,var,ell,eps,N,zeta)});
    f = f_map(const_id);
    phiextsol = fsolve(f, phiextinit, options);
    uold = [phiextsol; var];
    
    % calculate second point
    var=var+dvar;
    f_map = containers.Map(const_ids, {@(phiext) int_eq2d(phiext,var,const,ell,eps,N,zeta), @(phiext) int_eq2d(phiext,const,var,ell,eps,N,zeta)});
    f = f_map(const_id);
    phiextsol = fsolve(f, phiextsol, options);
    u0 = [phiextsol; var];

    % calculate first secant
    sec = u0 - uold;
    sec = sec/norm(sec);

    % secant continuation
    for j = 1:numiter
        uinit = u0 + sec*ds;
        fext = @(u) int_eq2d_ext(u,const_id,const,ell,eps,N,zeta,sec,uinit);   
        [usol, ~, ~, output] = fsolve(fext, uinit, options);
        if regrid_error >= 0 % check if regridding enabled
            err=abs(fft(usol(1:end-1)));
            err=max(err(3*N/8:5*N/8));
            if err > regrid_error % regrid
                disp('Now refining grid to new N = ' + string(N*2) + ' for iteration ' + string(j))
                
                usolh=fft(usol(1:end-2));
                usolh = [usolh(1:N/2);zeros(N,1);usolh(1+N/2:N)];
                usol = [2*ifft(usolh,'symmetric');usol(N+1:N+2)];
                
                u0h = fft(u0(1:end-2));
                u0h = [u0h(1:N/2);zeros(N,1);u0h(1+N/2:N)];
                u0 = [2*ifft(u0h,'symmetric');u0(N+1:N+2)];
                
                sec = usol - u0;
                uinit = u0 + sec*ds;
                
                N=2*N;
                ell = [0:N/2 -N/2+1:-1]'; % variable in Fourier space
                zeta = (0:L/N:L-L/N)';
                fext = @(u) int_eq2d_ext(u,const_id,const,ell,eps,N,zeta,sec,uinit);
                [usol, ~, ~, output] = fsolve(fext,usol, options);
            end
        end
        kx_array(2*i-1:2*i, j) = usol(end-1:end);
        sec = usol - u0;
        sec = sec/norm(sec);
        u0=usol;
        
        if adaptive_ds % adjust ds accordingly for next iteration
            if output.iterations > 8
                ds = ds/2;
                disp('ds decreased to ' + string(ds) + 'for iteration ' + string(j+1))
            elseif output.iterations < 2
                ds = ds*2;
                disp('ds increased to ' + string(ds) + 'for iteration ' + string(j+1))
            else
                % don't change ds
            end
        end
    end
end
end