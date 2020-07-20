function DF = jac_inteq2d(dphiext,phiext,var,dvar,const,ell,eps,N,zeta)
% for adding Jacobian to fsolve
% currently only for k_x vs. c_x with k_y fixed

phi = phiext(1:end-1); % function
k_x = phiext(end); % k_x variable as given in the notes
h=1+eps*sin(phi-zeta)-k_x; % g(u)-k_x in the notes

dphi = dphiext(1:end-1);
dk_x = dphiext(end);
% dh = eps*cos(phi-zeta).*dphi-dk_x;

% define 2D derivative operator in Fourier space
if var == 0
  Der = -abs(ell*const);
else
  Der = var/2-sqrt(var^2/4+(ell.^2)*const^2+k_x*var*1i*ell);
end

% F0 = phi+ifft((Der-1).^(-1).*fft(phi-h),'symmetric'); % from int_eq2d
% assuming var != 0, var = c_x
dphi_F0 = dphi + ifft((Der-1).^(-1).*fft(dphi - eps*cos(phi-zeta).*dphi), 'symmetric'); % (dF0/dphi)*dphi

Der_kx = (var*1i*ell)./((Der-1).^2.*2.*sqrt(var^2/4+(ell.^2)*const^2+k_x*var*1i*ell)); % d/dkx (Der-1).^(-1)
dkx_F0 = ifft(Der_kx.*fft(phi-h)*dk_x + (Der-1).^(-1).*fft(dk_x*ones(N, 1)), 'symmetric'); % (dF0/dk_x)*dk_x

Der_var = (1/2 - (var/2+k_x*1i*ell)./(2*sqrt(var^2/4+(ell.^2)*const^2+k_x*var*1i*ell)))./((Der-1).^2)*(-1); % d/dvar (Der-1).^(-1)
dvar_F0 = ifft(Der_var.*fft(phi-h)*dvar, 'symmetric'); % (dF0/dvar)*dvar

DF = [dphi_F0 + dkx_F0 + dvar_F0; sum(dphi)/N];



% 
% 
% 
% 
% dvar_Der = 1/2 - (var/2 + k_x*1i*ell)./(2*sqrt(var^2/4+(ell.^2)*const^2+k_x*var*1i*ell));
% dvar_Derinv = -dvar_Der./Der.^2;
% 
% DF = [dphi+ifft((Der-1).^(-1).*fft(dphi-dh),'symmetric'); sum(dphi)/N];
% Dvar= ifft(dvar_Derinv.*fft(phi-h)*dvar, 'symmetric');
% DF=DF+[Dvar; 0];
end