% main function where constants are defined and other functions are called

%close all;
clear all;

% define optimization parameters
N = 2^9; %  number of Fourier modes
numiter = 100; % number of secant continuation iterations
ds = 0.8; % secant continuation step size
dvar = 1e-4; % "baby continuation" step size/direction
options = optimset('Jacobian','off','Display','iter','TolFun',1e-6,'TolX',1e-6,'MaxIter',10,'Algorithm','trust-region-reflective');
regrid_error = 1e-4; % changes how fsolve_phi runs

% define constants
eps = 0.3; % theta parameter
c_x = 0.0;
%c_xs = 0.1; % fixed values for fsolve_kx only
k_y = 0.2;
%k_ys = 1; % fixed values for fsolve_kx only
L = 2*pi; % length of domain
k_xinit = 1; % initial k_x; notice that k_x = -c_y*k_y/c_x with c_y = -c_x/k_y
phiinit = zeros(N,1); % initial function
ell = [[0:N/2] [-N/2+1: -1]]'; % variable in Fourier space
zeta = (0:L/N:L-L/N)';



% NEW CODE
%{
k_yarray = [k_y:k_y:50*k_y]; % obtaining k_y values from the above continuation
sumarray = zeros(length(k_yarray),1);
for m = 1:length(k_yarray)
  [zetagrid, phigrid, fftphigrid, cgrid] = fsolve_phi(N, L, ell, zeta, eps, c_x, k_yarray(m), options, -1); % get fftphigrid from each value k_y
  sumarray(m) = sum(abs(ell).*(abs(fftphigrid).^2));
end
%}
% What I think it does: takes a grid of k_y values (could be anything,
% ideally greater than k_y = 0.2 since that has worked well in the past)
% and, for each k_y value, obtains a plot of the "y-values" of the Fourier 
% transform of phi (the periodic v in previous notes), and sums up those
% y-values with the corresponding entries in ell. Then, sumarray are all
% the sums which act as "y-values" for the plot being created.  The
% "x-values" to be plotted are simply the list of k_y's in k_yarray.



% ========== sample code for using all of the other functions ==========
% % plot u = phi - zeta at x = 0
% [zeta, phigrid, fftphigrid, kxgrid, cgrid] = fsolve_phi(N, L, ell, zeta, eps, c_x, k_y, options, -1);
% figure('Name', 'plot_u')
% plot_u(zeta, phigrid, cgrid);
% 
% % plot FT(phi) vs. zeta
% figure('Name', 'plot_ftphi')
% plot_ftphi(zeta, fftphigrid, cgrid);
% 
% % plot k_x vs. c_x for each k_y in k_ys
% [fsk1, fsk2, fsk3, kx_array] = fsolve_kx("k_y", k_ys, "c_x", c_x, N, L, numiter, ell, zeta, eps, k_xinit, phiinit, ds, dvar, options, -1, 0);
% figure('Name', 'plot_kx in c_x')
% plot_kx(fsk1, fsk2, fsk3, kx_array, [], [], [])
%
% % as above, with monotone plots and global minima of non-monotone plots marked
% monoinc = seqtype_kx(1, fsk4, fsk2);
% monodec = seqtype_kx(1i, fsk4, fsk2);
% [~, ~, pts] = seqtype_kx(0, fsk4, fsk2);
% figure('Name', 'plot_kx in c_x, monotone marked')
% plot_kx(fsk1, fsk2, fsk3, fsk4, monoinc, monodec, pts)
% 
% % for non-monotone plots of k_x vs. c_x, plot the (k_y, c_x) corresponding to global minima
% figure('Name', 'plot_extrloc')
% plot_extrloc(fsk2, fsk4, pts)
% 
% % plot k_x vs. k_y for each c_x in c_xs
% [fsk5, fsk6, fsk7, fsk8] = fsolve_kx("c_x", c_xs, "k_y", k_y, N, numiter, ell, zeta, eps, k_xinit, phiinit, ds, dvar, options);
% figure('Name', 'plot_kx in k_y')
% plot_kx(fsk5, fsk6, fsk7, fsk8, [], [], [])
% 
% % given k_y, plot k_x vs. c_x with regridding on
% [zeta, k_xgrid, ~, cgrid] = fsolve_phi(N, L, ell, zeta, eps, c_x, k_y, options, regrid_error);
% figure('Name', "fsolve_phi: k_x vs. c_x for k_y = " + string(k_y))
% plot_kx("k_y", [k_y], "c_x", [k_xgrid; cgrid], [], [], [])
% ======================================================================


% given k_y, plot k_x vs. c_x with regridding on
% [zeta, k_xgrid, ~, cgrid] = fsolve_phi(N, L, ell, zeta, eps, c_x, k_y, options, regrid_error);
% figure('Name', "fsolve_phi: k_x vs. c_x for k_y = " + string(k_y))
% plot_kx("k_y", [k_y], "c_x", [k_xgrid; cgrid], [], [], [])

% % plot k_x vs. c_x for each k_y in k_ys
% % dotted line = monotone increasing, dot-dashed line = monotone decreasing
% [const_id, constlist, var_id, kx_array] = fsolve_kx("k_y", k_ys, "c_x", c_x, N, numiter, ell, zeta, eps, k_xinit, phiinit, ds, dvar, options);
% monoinc = seqtype_kx(1, kx_array, constlist);
% monodec = seqtype_kx(1i, kx_array, constlist);
% [~, ~, extr_pts] = seqtype_kx(0, kx_array, constlist);
% figure('Name', 'plot_kx in c_x')
% plot_kx(const_id, constlist, var_id, kx_array, monoinc, monodec, extr_pts)
% 
% % for non-monotone plots of k_x vs. c_x, plot the (k_y, c_x) corresponding to global minima
% figure('Name', 'plot_extrloc')
% plot_extrloc(constlist, kx_array, extr_pts)

[zeta, phigrid, fftphigrid, kxgrid, cgrid] = fsolve_phi(N, L, ell, zeta, eps, c_x, k_y, options, -1);
Der = (cgrid./2-sqrt((cgrid.^2)./4+(ell.^2).*k_y^2+kxgrid.*cgrid.*1i.*ell));
colors = distinguishable_colors(length(cgrid));

vid = VideoWriter(['PhaseDiffusionSolution'],'MPEG-4');
vid.Quality = 100;
Framerate = 120;
open(vid);

%{
% Plotting u(x,zeta) for each x in an animation
figure('Name', 'plot_ufull2');
for x = 0:0.1:20
    sol = ifft(exp(Der.*x).*fft(phigrid),'symmetric');
    plot(0);
    for m = 1:length(cgrid)
        hold on
        plot(zeta,sol(:,m)-zeta+kxgrid(m)*x, 'Color', colors(m, :))
        hold off
    end
    
    %title(legend(string(cgrid)), 'c')
    %legend("Location", "eastoutside")
    xlabel('zeta')
    ylabel('u(x,zeta)')
    title(['Phase-Diffusion, x = ' num2str(x)])
    axis([0 2*pi -2*pi 20])
    
    drawnow;
    frame = getframe(gcf);
    writeVideo(vid,frame);
end
%}


% Plotting 2D figure u(x,zeta)
figure('Name', 'plot_ufull');
set(gcf, 'position', [100,100,1000,9*6.3*10])
matrix = zeros(N,100);
for x = 0:1:length(matrix(1,:))-1
    sol = ifft(exp(Der.*x./10).*fft(phigrid),'symmetric')-zeta+kxgrid*(x/10);
    matrix(:,x+1) = sol;
end
set(gca,'YDir','reverse')
imagesc(flipud(matrix)); colormap parula; % initially set to "bone"; set to "jet" for Montie's pictures
colorbar;
xlabel('x')
ylabel('zeta')
    
A = [10:10:100];
B = [1:63.8:512];
xticks(A);
yticks(B);
Xticks = {};
Yticks = {};
for i = 1:length(A)
    Xticks{i} = round((length(matrix(1,:))/10)*i/length(A),1);
end
xticklabels(Xticks);
for i = 1:length(B)
    %Yticks{length(B)-i+1} = round(2*pi*i/length(B),1);
    Yticks{i} = round((2*pi/8)*(9-i),1);
end
yticklabels(Yticks);

%drawnow;
%frame = getframe(gcf);
%writeVideo(vid,frame);
%}

close(vid);
%}




