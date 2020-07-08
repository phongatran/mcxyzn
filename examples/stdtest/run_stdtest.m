clear
addpath ../..
addpath ../../matlab_functions/

%%% Basic MC configuration %%%
cfg.SAVEON = 1; % 1 = save myname_T.bin, myname_H.mci 
                % 0 = don't save. Just check the program.
cfg.name = 'stdtestM';
cfg.time = 1.0;               %Simulation time in min
Nbins = 200;
cfg.binsize = 10/Nbins;        %Length of a voxel
cfg.dim = [Nbins,Nbins,Nbins]; %Number of voxels in each direction [Nx,Ny,Nz]

cfg.mcflag       = 0;   % launch: 0 = uniform beam, 1 = Gaussian, 2 = isotropic pt. 
                        % 3 = rectangular beam (use xfocus,yfocus for x,y halfwidths)
cfg.launchflag   = 0;   % 0 = let mcxyz.c calculate launch trajectory
                        % 1 = manually set launch vector.
cfg.boundaryflag = 1;   % 0 = no boundaries, 1 = escape at boundaries
                        % 2 = escape at surface only. No x,y, bottom z
                        % boundaries
cfg.gradientflag = 1;   % 0 = Fresnel's law using the voxel faces [inaccurate for curved and oblique geometries]
                        % 1 = gradient-based Fresnel's laws using the Sobel
                        % filter
                        % 2 = gradient-based Fresnel's laws + additional smoothing filter 
cfg.srcpos   = [0,0,0];   % Set position of source [x,y,z];
cfg.srcfocus = [0,0,inf]; % Set position of focus, so mcxyz can calculate launch trajectory
                          % [xfocus,yfocus,zfocus];
cfg.radius   = 0;         % 1/e radius of beam at tissue surface
cfg.waist    = 0;         % 1/e radius of beam at focus

% only used if launchflag == 1 (manually set launch trajectory) / ux^2 + uy^2 + uz^2 = 1
cfg.launchvec = [0.7,0.4,sqrt(1 - 0.7^2 - 0.4^2)]; 

cfg.Nt = 2;
i=1; % matched non-scattering medium
cfg.muav(i)  = 0.001;
cfg.musv(i)  = 100;
cfg.gv(i)    = 1.0; % <---- non-scattering
cfg.nv(i)    = 1.4;
i=2; % standard tissue
cfg.muav(i)  = 1;
cfg.musv(i)  = 100;
cfg.gv(i)    = 0.90;
cfg.nv(i)    = 1.4;

%%% Setting up simulation volume %%% 
T = double(zeros(cfg.dim(1),cfg.dim(2),cfg.dim(3))); 
T = T + 2;      % fill background w standard tissue
zsurf = 0.500;  % position of air/skin surface
izsurf = round(zsurf/cfg.binsize);
T(:,:,1:izsurf) = 1; % top layer nonscattering

cfg.T = T;

create_simfiles(cfg);
launch_simulation(cfg);
%[fluence] = plot_mcresults(cfg);


