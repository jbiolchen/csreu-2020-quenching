function F = int_eq2d(phiext,c_x,k_y,ell,eps,N,zeta)
% applies boundary integral operator to phiext

phi = phiext(1:end-1); % function
k_x = phiext(end); % k_x variable as given in the notes
h = 1+eps*sin(phi-zeta)-k_x; % g(u)-k_x in the notes

% define 2D derivative operator in Fourier space
Der = c_x/2-sqrt(c_x^2/4+(ell.^2)*k_y^2+k_x*c_x*1i*ell);

F = [phi+ifft((Der-1).^(-1).*fft(phi-h),'symmetric'); sum(phi)/N];
end
