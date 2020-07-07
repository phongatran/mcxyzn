clear
addpath ../..
addpath ../../matlab_functions/
[cfg,F,T] = readfluence('lens');


SAVEPICSON = 0;
if SAVEPICSON
    sz = 10; fz = 7; fz2 = 5; % to use savepic.m
else
    sz = 12; fz = 9; fz2 = 7; % for screen display
end

Nx = cfg.dim(1); Ny = cfg.dim(2); Nz = cfg.dim(3);
xs = cfg.srcpos(1); ys = cfg.srcpos(2); zs = cfg.srcpos(3);
xfocus = cfg.srcfocus(1); yfocus = cfg.srcfocus(2); zfocus = cfg.srcfocus(3);
nm = 660;

x = ([1:Nx]-Nx/2-1/2)*cfg.binsize;
y = ([1:Ny]-Ny/2-1/2)*cfg.binsize;
z = ([1:Nz]-1/2)*cfg.binsize;
ux = [2:Nx-1];
uy = [2:Ny-1];
uz = [2:Nz-1];
zmin = min(z);
zmax = max(z);
zdiff = zmax-zmin;
xmin = min(x);
xmax = max(x);
xdiff = xmax-xmin;

Tzx = reshape(T(:,Ny/2,:),Nx,Nz)';
PRINTON = 0;
tissue = makeTissueList(nm);
Nt = length(tissue);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot tissue structure %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);clf
imagesc(x(ux),z(uz),Tzx(uz,ux))
hold on
cmap = makecmap(Nt);
colormap(cmap)
colorbar
set(gca,'fontsize',sz)
set(colorbar,'fontsize',1)
xlabel('x [cm]')
ylabel('z [cm]')
% title('Tissue','fontweight','normal','fontsize',fz2)
% for i=1:Nt
%     yy = zmin + (Nt-i)/(Nt-1)*zdiff;
%     text(xmax*1.4,yy, sprintf('%d %s',i,tissue(i).name),'fontsize',fz2)
% end

% draw launch
N = 20; % # of beam rays drawn
switch cfg.mcflag
    case 0 % uniform
        for i=0:N
            plot((-cfg.radius + 2*cfg.radius*i/N)*[1 1],[zs max(z)],'r-')
        end
    case 1 % Gaussian
        for i=0:N
            plot([(-cfg.radius + 2*cfg.radius*i/N) xfocus],[zs zfocus],'r-')
        end
    case 2 % iso-point
        for i=1:N
            th = (i-1)/19*2*pi;
            xx = Nx/2*cos(th) + xs;
            zz = Nx/2*sin(th) + zs;
            plot([xs xx],[zs zz],'r-')
        end
    case 3 % rectangle
        zz = max(z);
        for i=0:N
            xx = -radius + 2*radius*i/N;
            plot([xx xx],[zs zz],'r-')
        end
end
axis equal image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot fluence in the zx plane %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fzx = reshape(F(:,Ny/2,:),Nx,Nz)'; % in z,x plane through source

figure(2);clf
imagesc(x,z,log10(Fzx))
hold on
text(max(x)*1.2,min(z)-0.04*max(z),'log_{10}( \phi )','fontsize',fz)
colorbar
set(gca,'fontsize',sz)
xlabel('x [cm]')
ylabel('z [cm]')
title('Fluence \phi [W/cm^2/W.delivered] ','fontweight','normal','fontsize',fz)
colormap(makec2f)
axis equal image
%axis([min(x) max(x) min(z) max(z)])
text(min(x)-0.2*max(x),min(z)-0.08*max(z),sprintf('runtime = %0.1f min',cfg.time),...
    'fontsize',fz2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot fluence in the zx plane %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fzy = reshape(F(Nx/2,:,:),Ny,Nz)'; % in z,x plane through source

figure(3);clf
imagesc(x,z,log10(Fzy))
hold on
text(max(x)*1.2,min(z)-0.04*max(z),'log_{10}( \phi )','fontsize',fz)
colorbar
set(gca,'fontsize',sz)
xlabel('y [cm]')
ylabel('z [cm]')
title('Fluence \phi [W/cm^2/W.delivered] ','fontweight','normal','fontsize',fz)
colormap(makec2f)
axis equal image
%axis([min(x) max(x) min(z) max(z)])
text(min(x)-0.2*max(x),min(z)-0.08*max(z),sprintf('runtime = %0.1f min',cfg.time),...
    'fontsize',fz2)
