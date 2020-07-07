function [cfg,F,T] = readfluence(name)
% lookmcxyz.m
%   if quickON==1, no figures. Must specify myname.
%
%   Looks at myname_F.bin, created by mcxyz.c 
%   where myname is the name of the run: myname_T.bin, myname_H.mci
%   Makes figures:
%       myname_tissue.jpg   = tissue structure (shows tissue types)
%       myname_Fzx.jpg      = fluence rate vs z,x
%       myname_Fzy.jpg      = fluence rate vs z,y
%   Uses:
%       myname_H.mci    = input file from maketissue.m
%       myname_T.bin    = tissue input file from maketissue.m
%       myname_F.bin    = fluence rate output from Monte Carlo
%       reportH_mci.m   = lists input parameters in myname_H.mci
%       makecmap.m      = makes colormap for tissue types
%       makec2f.m       = makes colormap for fluence rate
%
%   This example sets myname = 'skinvessel'.
%
% 7/feb/2017, add boundaryflag (see A(10)).
% 1/june/2017 , no major changes, just clean up display outputs.
% Steven L Jacques
%home; clear

% Load header file
filename = sprintf('%s_H.mci',name);
disp(['loading ' filename])
fid = fopen(filename, 'r');
A = fscanf(fid,'%f',[1 Inf])';
fclose(fid);

%% parameters
cfg.time = A(1);
cfg.dim(1) = A(2);
cfg.dim(2) = A(3);
cfg.dim(3) = A(4);
cfg.binsize = A(5);
%dy = A(6);
%dz = A(7);
cfg.mcflag = A(8);
cfg.launchflag = A(9);
cfg.boundaryflag = A(10);
cfg.gradientflag = A(11);
cfg.srcpos(1) = A(12);
cfg.srcpos(2) = A(13);
cfg.srcpos(3) = A(14);
cfg.srcfocus(1) = A(15);
cfg.srcfocus(2) = A(16);
cfg.srcfocus(3) = A(17);
cfg.launchvec(1) = A(18);
cfg.launchvec(2) = A(19);
cfg.launchvec(3) = A(20);
cfg.radius = A(21);
cfg.waist = A(22);
cfg.Nt = A(23);
j = 23;
for i=1:cfg.Nt
    j=j+1;
    cfg.muav(i,1) = A(j);
    j=j+1;
    cfg.musv(i,1) = A(j);
    j=j+1;
    cfg.gv(i,1) = A(j);
    j=j+1;
    cfg.nv(i,1) = A(j);
end

reportHmci(name);


Nx = cfg.dim(1);
Ny = cfg.dim(2);
Nz = cfg.dim(3);

%% Load Fluence rate F(x,y,z) 
filename = sprintf('%s_F.bin',name);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Nx*Ny*Nz, 'float');
    fclose(fid);
toc
F = reshape(Data,Nx,Ny,Nz); % F(x,y,z)
%F = permute(F,[2,1,3]);
%%
% Load tissue structure in voxels, T(x,y,z) 
filename = sprintf('%s_T.bin',name);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Nx*Ny*Nz, 'uint8');
    fclose(fid);
toc
T = reshape(Data,Nx,Ny,Nz); % T(x,y,z)
%T = permute(T,[2,1,3]);
end

%%
% x = ([1:Nx]-Nx/2-1/2)*dx;
% y = ([1:Ny]-Ny/2-1/2)*dx;
% z = ([1:Nz]-1/2)*dz;
% ux = [2:Nx-1];
% uy = [2:Ny-1];
% uz = [2:Nz-1];
% zmin = min(z);
% zmax = max(z);
% zdiff = zmax-zmin;
% xmin = min(x);
% xmax = max(x);
% xdiff = xmax-xmin;

% format compact
% commandwindow
% 
% quickON = 0
% 
% if ~quickON
%     %%%% USER CHOICES <---------- you must specify -----
%     %
%     nm = 660;
%     %%%%
% 
%     SAVEPICSON = 0;
%     if SAVEPICSON
%         sz = 10; fz = 7; fz2 = 5; % to use savepic.m
%     else
%         sz = 12; fz = 9; fz2 = 7; % for screen display
%     end
%     disp(sprintf('------ mcxyz %s -------',name))
% end

%% Look at structure, Tzx

% 
% %% Look at Fluence Fzx @ launch point
% Fzx = reshape(F(Ny/2,:,:),Nx,Nz)'; % in z,x plane through source
% 
% figure(2);clf
% imagesc(x,z,log10(Fzx))
% hold on
% text(max(x)*1.2,min(z)-0.04*max(z),'log_{10}( \phi )','fontsize',fz)
% colorbar
% set(gca,'fontsize',sz)
% xlabel('x [cm]')
% ylabel('z [cm]')
% title('Fluence \phi [W/cm^2/W.delivered] ','fontweight','normal','fontsize',fz)
% colormap(makec2f)
% axis equal image
% %axis([min(x) max(x) min(z) max(z)])
% text(min(x)-0.2*max(x),min(z)-0.08*max(z),sprintf('runtime = %0.1f min',time_min),...
%     'fontsize',fz2)
% 
% %% look Fzy
% Fzy = reshape(F(:,Nx/2,:),Ny,Nz)';
% 
% iy = round((dy*Ny/2 + 0.15)/dy);
% iz = round(zs/dz);
% zzs  = zs;
% %Fdet = mean(reshape(Fzy(iz+[-1:1],iy+[0 1]),6,1));
% 
% figure(3);clf
% imagesc(y,z,log10(Fzy))
% hold on
% text(max(x)*1.2,min(z)-0.04*max(z),'log_{10}( \phi )','fontsize',fz)
% colorbar
% set(gca,'fontsize',sz)
% xlabel('y [cm]')
% ylabel('z [cm]')
% title('Fluence \phi [W/cm^2/W.delivered] ','fontweight','normal','fontsize',fz)
% colormap(makec2f)
% axis equal image
% text(min(x)-0.2*max(x),min(z)-0.08*max(z),sprintf('runtime = %0.1f min',time_min),...
%     'fontsize',fz2)
% 
% 
% %% look Azx
% Fzx = reshape(F(Ny/2,:,:),Nx,Nz)'; % in z,x plane through source
% mua = muav(reshape(T(Ny/2,:,:),Nx,Nz)');
% Azx = Fzx.*mua;
% 
% figure(2);clf
% imagesc(x,z,log10(Azx))
% hold on
% text(max(x)*1.2,min(z)-0.04*max(z),'log_{10}( A )','fontsize',fz)
% colorbar
% set(gca,'fontsize',sz)
% xlabel('x [cm]')
% ylabel('z [cm]')
% title('Deposition A [W/cm^3/W.delivered] ','fontweight','normal','fontsize',fz)
% colormap(makec2f)
% axis equal image
% %axis([min(x) max(x) min(z) max(z)])
% text(min(x)-0.2*max(x),min(z)-0.08*max(z),sprintf('runtime = %0.1f min',time_min),...
%     'fontsize',fz2)
% 
% drawnow
