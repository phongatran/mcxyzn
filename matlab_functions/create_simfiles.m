function [] = create_simfiles(cfg)

% maketissue_stdtestM.m
% maketissue.m
%   Creates a cube of optical property pointers,T(y,x,z), saved in
%       myname_T.bin = a tissue structure file
%   which specifies a complex tissue for use by mcxyz.c.
%
%   Also prepares a listing of the optical properties at chosen wavelength
%   for use by mcxyz.c, [mua, mus, g], for each tissue type specified
%   in myname_T.bin. This listing is saved in
%       myname_H.mci = the input file for use by mcxyz.c.
%
%   Will generate a figure illustrating the tissue with its various
%   tissue types and the beam being launched.
%
%   Uses
%       makeTissueList.m
%
%   To use, 
%       1. Prepare makeTissueList.m so that it contains the tissue
%   types desired.
%       2. Specify the USER CHOICES.
%       2. Run this program, maketissue.m.
%
%   Note: mcxyz.c can use optical properties in cm^-1 or mm^-1 or m^-1,
%       if the bin size (binsize) is specified in cm or mm or m,
%       respectively.
%
%  Steven L. Jacques. updated Aug 21, 2014.
%       

% clear
% format compact
% clc
% home
         
% Sets position of source
xs          = cfg.srcpos(1);      	% x of source
ys          = cfg.srcpos(2);        % y of source
zs          = cfg.srcpos(3);  	% z of source

% Set position of focus, so mcxyz can calculate launch trajectory
xfocus      = cfg.srcfocus(1);        % set x,position of focus
yfocus      = cfg.srcfocus(2);        % set y,position of focus
zfocus      = cfg.srcfocus(3);    	% set z,position of focus (=inf for collimated beam)

% only used if mcflag == 0 or 1 or 3 (not 2=isotropic pt.)
radius      = cfg.radius;   % 1/e radius of beam at tissue surface
waist       = cfg.waist;  	% 1/e radius of beam at focus

% only used if launchflag == 1 (manually set launch trajectory):
ux0         = cfg.launchvec(1);      % trajectory projected onto x axis
uy0         = cfg.launchvec(2);      % trajectory projected onto y axis
uz0         = cfg.launchvec(3); % such that ux^2 + uy^2 + uz^2 = 1
%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% 
% Prepare Monte Carlo 
%%%%%%%%%%


    
% Specify Monte Carlo parameters    
Nx = cfg.dim(1);
Ny = cfg.dim(2);
Nz = cfg.dim(3);
dx = cfg.binsize;
dy = cfg.binsize;
dz = cfg.binsize;
x  = ([1:Nx]'-Nx/2)*dx;
y  = ([1:Ny]'-Ny/2)*dy;
z  = [1:Nz]'*dz;
zmin = min(z);
zmax = max(z);
xmin = min(x);
xmax = max(x);

if isinf(zfocus), zfocus = 1e12; end

%%%%%%
% CREATE TISSUE STRUCTURE T(y,x,z)
%   Create T(y,x,z) by specifying a tissue type (an integer)
%   for each voxel in T.
%
%   Note: one need not use every tissue type in the tissue list.
%   The tissue list is a library of possible tissue types.

T = cfg.T;

%%
if cfg.gradientflag > 0
    unique_nv = unique(cfg.nv);
    gradient_map_x = zeros(size(T,1),size(T,2),size(T,3));
    gradient_map_y = zeros(size(T,1),size(T,2),size(T,3));
    gradient_map_z = zeros(size(T,1),size(T,2),size(T,3));

    for j = 1:length(unique_nv)
        n_map = zeros(size(T,1),size(T,2),size(T,3));
        for i = 1:length(n_map(:))
            if cfg.nv(T(i)) == cfg.nv(j)
                n_map(i) = 1;
            end
        end
        [Gx,Gy,Gz] = imgradientxyz(n_map);
        smooth_map = smoothn({Gx,Gy,Gz},2);
        gradient_map_x(n_map == 1) = smooth_map{1}(n_map==1);
        gradient_map_y(n_map == 1) = smooth_map{2}(n_map==1);
        gradient_map_z(n_map == 1) = smooth_map{3}(n_map==1);
        magn = sqrt(gradient_map_x.^2+gradient_map_y.^2+gradient_map_z.^2);
        gradient_map_x(magn < 1) = 0;
        gradient_map_y(magn < 1) = 0;
        gradient_map_z(magn < 1) = 0;
    end
    gradient_magn = sqrt(gradient_map_x.^2+gradient_map_y.^2+gradient_map_z.^2);
    gradient_map_x = gradient_map_x./gradient_magn;
    gradient_map_y = gradient_map_y./gradient_magn;
    gradient_map_z = gradient_map_z./gradient_magn;
end

if cfg.SAVEON
    tic
    % convert T to linear array of integer values, v(i)i = 0;
    v = uint8(reshape(T,Ny*Nx*Nz,1));
 
    %% WRITE FILES
    % Write myname_H.mci file
    %   which contains the Monte Carlo simulation parameters
    %   and specifies the tissue optical properties for each tissue type.
    commandwindow
    disp(sprintf('--------create %s --------',cfg.name))
    filename = sprintf('%s_H.mci',cfg.name);
    fid = fopen(filename,'w');
        % run parameters
        fprintf(fid,'%0.2f\n',cfg.time);
        fprintf(fid,'%d\n'   ,Nx);
        fprintf(fid,'%d\n'   ,Ny);
        fprintf(fid,'%d\n'   ,Nz);
        fprintf(fid,'%0.4f\n',dx);
        fprintf(fid,'%0.4f\n',dy);
        fprintf(fid,'%0.4f\n',dz);
        % launch parameters
        fprintf(fid,'%d\n'   ,cfg.mcflag);
        fprintf(fid,'%d\n'   ,cfg.launchflag);
        fprintf(fid,'%d\n'   ,cfg.boundaryflag);
        fprintf(fid,'%d\n'   ,cfg.gradientflag);
        fprintf(fid,'%0.4f\n',xs);
        fprintf(fid,'%0.4f\n',ys);
        fprintf(fid,'%0.4f\n',zs);
        fprintf(fid,'%0.4f\n',xfocus);
        fprintf(fid,'%0.4f\n',yfocus);
        fprintf(fid,'%0.4f\n',zfocus);
        fprintf(fid,'%0.4f\n',ux0); % if manually setting ux,uy,uz
        fprintf(fid,'%0.4f\n',uy0);
        fprintf(fid,'%0.4f\n',uz0);
        fprintf(fid,'%0.4f\n',radius);
        fprintf(fid,'%0.4f\n',waist);
        % tissue optical properties
        fprintf(fid,'%d\n',cfg.Nt);
        for i=1:cfg.Nt
            fprintf(fid,'%0.4f\n',cfg.muav(i));
            fprintf(fid,'%0.4f\n',cfg.musv(i));
            fprintf(fid,'%0.4f\n',cfg.gv(i));
            fprintf(fid,'%0.4f\n',cfg.nv(i));
        end
    fclose(fid);
 
    %% write myname_T.bin file
    filename = sprintf('%s_T.bin',cfg.name);
    disp(['create ' filename])
    fid = fopen(filename,'wb');
    fwrite(fid,v,'uint8');
    fclose(fid);
    
    if cfg.gradientflag > 0
        %% write myname_Gx.bin file
        gradient_map_x = reshape(gradient_map_x,Ny*Nx*Nz,1);
        filename = sprintf('%s_Gx.bin',cfg.name);
        disp(['create ' filename])
        fid = fopen(filename,'wb');
        fwrite(fid,gradient_map_x,'float');
        fclose(fid);
        %% write myname_Gy.bin file
        gradient_map_y = reshape(gradient_map_y,Ny*Nx*Nz,1);
        filename = sprintf('%s_Gy.bin',cfg.name);
        disp(['create ' filename])
        fid = fopen(filename,'wb');
        fwrite(fid,gradient_map_y,'float');
        fclose(fid);

        %% write myname_Gz.bin file
        gradient_map_z = reshape(gradient_map_z,Ny*Nx*Nz,1);
        filename = sprintf('%s_Gz.bin',cfg.name);
        disp(['create ' filename])
        fid = fopen(filename,'wb');
        fwrite(fid,gradient_map_z,'float');
        fclose(fid);
        toc
    end
end % SAVEON
 

%% Look at structure of Tzx at iy=Ny/2
%Txzy = shiftdim(T,1);   % Tyxz --> Txzy
Txz = reshape(T(:,Ny/2,:),[Nx Nz]);
Tzx = shiftdim(Txz,1);

%%
figure(1); clf
sz = 12;  fz = 10; 
imagesc(x,z,Tzx,[1 cfg.Nt])
hold on
set(gca,'fontsize',sz)
xlabel('x [cm]')
ylabel('z [cm]')
colorbar
cmap = makecmap(cfg.Nt);
colormap(cmap)
set(colorbar,'fontsize',1)
% label colorbar
zdiff = zmax-zmin;
%%%

% for i=1:Nt
%     yy = (Nt-i)/(Nt-1)*Nz*dz;
%     text(max(x)*1.2,yy, tissue(i).name,'fontsize',fz)
% end
% 
text(xmax,zmin - zdiff*0.06, 'Tissue types','fontsize',fz)
axis equal image
axis([xmin xmax zmin zmax])

%%% draw launch
N = 20; % # of beam rays drawn
switch cfg.mcflag
    case 0 % uniform
        for i=0:N
            plot((-radius + 2*radius*i/N)*[1 1],[zs max(z)],'r-')
        end

    case 1 % Gaussian
        for i=0:N
            plot([(-radius + 2*radius*i/N) xfocus],[zs zfocus],'r-')
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
        for i=1:N
            xx = -radius + 2*radius*i/20;
            plot([xx xx],[zs zz],'r-')
        end
end

disp('done')
end
