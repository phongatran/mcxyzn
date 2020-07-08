clear
addpath ../..
addpath ../../matlab_functions/

%%% Basic MC configuration %%%
cfg.SAVEON = 1; % 1 = save myname_T.bin, myname_H.mci 
                % 0 = don't save. Just check the program.
cfg.name = 'lens_water';
cfg.time = 1;               %Simulation time in min
Nbins = 200;
cfg.binsize = 0.01;        %Length of a voxel
cfg.dim = [Nbins,Nbins,Nbins]; %Number of voxels in each direction [Nx,Ny,Nz]

cfg.mcflag       = 0;   % launch: 0 = uniform beam, 1 = Gaussian, 2 = isotropic pt. 
                        % 3 = rectangular beam (use xfocus,yfocus for x,y halfwidths)
cfg.launchflag   = 0;   % 0 = let mcxyz.c calculate launch trajectory
                        % 1 = manually set launch vector.
cfg.boundaryflag = 1;   % 0 = no boundaries, 1 = escape at boundaries
                        % 2 = escape at surface only. No x,y, bottom z
                        % boundaries
cfg.gradientflag = 1;   % 0 = Fresnel's law using the voxel faces [inaccurate for curved and oblique geometries]
                        % 1 = gradient-based Fresnel's laws w. smoothing
                        % 2 = gradient-based Fresnel's laws w. smoothing and interpolation 
cfg.srcpos   = [0,0,0.0001];   % Set position of source [x,y,z];
cfg.srcfocus = [0,0,inf]; % Set position of focus, so mcxyz can calculate launch trajectory
                          % [xfocus,yfocus,zfocus];
cfg.radius   = 0.15;         % 1/e radius of beam at tissue surface
cfg.waist    = 0.15;         % 1/e radius of beam at focus

% only used if launchflag == 1 (manually set launch trajectory) / ux^2 + uy^2 + uz^2 = 1
cfg.launchvec = [0.7,0.4,sqrt(1 - 0.7^2 - 0.4^2)]; 

cfg.Nt   = 3;
cfg.muav = [0.0004,0.001,0.001];
cfg.musv = [1.0 1.0 100.0];
cfg.gv   = [1.0 1.0 1.0];
cfg.nv   = [1.52 1.33 1.33];

%%% Setting up simulation volume %%% 
T = double(zeros(cfg.dim(1),cfg.dim(2),cfg.dim(3))); 
zsurf = 0.0100;  % position of air/skin surface
for iz=1:cfg.dim(3) % for every depth z(iz)
    radius_c = 40/2;
    center = [200,200,170]/2;
    for i = 1:Nbins
        for j = 1:Nbins
            if (round(sqrt((center(1)-i)^2+(center(2)-j)^2+(center(3)-iz)^2))<radius_c)
               T(i,j,iz) = 1; 
            end
        end
    end
    radius_c = 50/2;
    center = [200,200,210]/2;
    for i = 1:Nbins
        for j = 1:Nbins
            if (round(sqrt((center(1)-i)^2+(center(2)-j)^2+(center(3)-iz)^2)<radius_c) && (T(i,j,iz)==1))
               T(i,j,iz) = 1; 
            elseif (iz>178/2)
               T(i,j,iz) = 3;
            else
               T(i,j,iz) = 2;
            end
        end
    end
end % iz
 
cfg.T = T;

create_simfiles(cfg);
launch_simulation(cfg);
%[fluence] = plot_mcresults(cfg);


