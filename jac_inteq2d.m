function DF = jac_inteq2d(dphiext,phiext,var,dvar,const_id,const,ell,eps,N,zeta)
% for adding Jacobian to fsolve

phi = phiext(1:end-1); % function
k_x = phiext(end); % k_x variable as given in the notes
h=1+eps*sin(phi-zeta)-k_x; % g(u)-k_x in the notes

dphi = dphiext(1:end-1);
dk_x = dphiext(end);

if const_id == "k_y"
    c_x = var;
    k_y = const;
else % const_id == "c_x"
    c_x = const;
    k_y = var;
end

% define 2D derivative operator in Fourier space
if c_x == 0
  Der = -abs(ell*k_y);
else
  Der = c_x/2-sqrt(c_x^2/4+(ell.^2)*k_y^2+k_x*c_x*1i*ell);
end

% define Jacobian info
% from int_eq2d, F0 = phi+ifft((Der-1).^(-1).*fft(phi-h),'symmetric');
dphi_F0 = dphi + ifft((Der-1).^(-1).*fft(dphi - eps*cos(phi-zeta).*dphi), 'symmetric'); % (dF0/dphi)*dphi
Der_kx = (c_x*1i*ell)./((Der-1).^2.*2.*sqrt(c_x^2/4+(ell.^2)*k_y^2+k_x*c_x*1i*ell)); % d/dkx (Der-1).^(-1)
if const_id == "k_y"
    Der_var = (1/2 - (c_x/2+k_x*1i*ell)./(2*sqrt(c_x^2/4+(ell.^2)*k_y^2+k_x*c_x*1i*ell)))./((Der-1).^2)*(-1); % d/dc_x (Der-1).^(-1)
else % const_id == "c_x"
    Der_var = ell.^2*k_y./(sqrt(c_x^2/4+(ell.^2)*k_y^2+k_x*c_x*1i*ell).*(Der-1).^2); % d/dk_y (Der-1).^(-1)
    if c_x == 0 % same as above, but NaN (i.e. 0/0) values replaced with 0
        Der_kx = 0;
        Der_var = abs(ell)./(Der-1).^2; % d/dk_y (Der-1).^(-1)
    end
end
dkx_F0 = ifft(Der_kx.*fft(phi-h)*dk_x + (Der-1).^(-1).*fft(dk_x*ones(N, 1)), 'symmetric'); % (dF0/dk_x)*dk_x
dvar_F0 = ifft(Der_var.*fft(phi-h)*dvar, 'symmetric'); % (dF0/dc_x)*dc_x or (dF0/dk_y)*dk_y

DF = [dphi_F0 + dkx_F0 + dvar_F0; sum(dphi)/N];
end