function [const_id, constlist, var_id, kx_array] = fsolve_kx(const_id, constlist, var_id, var, N, L, numiter, ell, zeta, eps, k_xinit, phiinit, ds, dvar, options, regrid_error, adaptive_ds)
% const is k_y or c_x; var is the other one
% for each const in constlist, calculate k_x as a function of var using secant continuation

% store input grid size to reset after each continuation
N_in = N;
ell_in = ell;
zeta_in = zeta;

const_ids = ["k_y", "c_x"]; % possible values of const_id
phiextinit = [phiinit;k_xinit]; % initial data

kx_array = zeros(2*length(constlist), numiter);
for i = 1:length(constlist)
    const = constlist(i);
    
    % reset to input grid size
    N = N_in;
    ell = ell_in;
    zeta = zeta_in;
    
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
        [usol, ~, flag, output] = fsolve(fext, uinit, options);
        
        % if adaptive ds enabled, do fsolve with smaller ds until it converges
        while adaptive_ds
            if flag > 0 % fsolve converged
                break
            end
            ds = ds/4;
            disp('Repeating iteration ' + string(j) + ' with smaller ds = ' + string(ds))
            uinit = u0 + sec*ds;
            fext = @(u) int_eq2d_ext(u,const_id,const,ell,eps,N,zeta,sec,uinit);   
            [usol, ~, flag, output] = fsolve(fext, uinit, options);
        end
        
        % if regridding enabled, do fsolve with larger N until err < regrid_error
        while regrid_error >= 0
            err=abs(fft(usol(1:end-2)));
            err=max(err(3*N/8:5*N/8));
            if err < regrid_error
                % using a stricter err definition, decrease N for iteration j+1 if reasonable
                err=abs(fft(usol(1:end-2)));
                err=max(err(1*N/8:7*N/8));
                if err < regrid_error
                    disp('Now decreasing grid size to new N = ' + string(N/2) + ' for iteration ' + string(j+1))
                    
                    usolh = fft(usol(1:end-2));
                    usolh = [usolh(1:N/4); usolh(3*N/4+1:N)];
                    usol = [ifft(usolh, 'symmetric')/2; usol(N+1:N+2)];
                    
                    u0h = fft(u0(1:end-2));
                    u0h = [u0h(1:N/4); u0h(3*N/4+1:N)];
                    u0 = [ifft(u0h,'symmetric')/2;u0(N+1:N+2)];

                    N=N/2;
                    ell = [0:N/2 -N/2+1:-1]'; % variable in Fourier space
                    zeta = (0:L/N:L-L/N)';
                end
                break
            end
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
        
        % if adaptive ds enabled, adjust ds to optimize fsolve for next iteration
        if adaptive_ds
            if output.iterations > 5
                ds = ds/2;
                disp('ds decreased to ' + string(ds) + 'for iteration ' + string(j+1))
            elseif output.iterations < 2
                ds = ds*1.3;
                disp('ds increased to ' + string(ds) + 'for iteration ' + string(j+1))
            else
                % don't change ds
            end
        end
        
        % store k_x, var for iteration j; update sec, u0 for iteration j+1
        kx_array(2*i-1:2*i, j) = usol(end-1:end);
        sec = usol - u0;
        sec = sec/norm(sec);
        u0=usol;  
    end
end
end