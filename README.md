# mcxyzn 
## Overview
*mcxyzn* is an extension of *mcxyz.c* with accurate handling of curved, complex, and oblique surfaces when accounting for refraction and reflection event when light hits a boundary with mismatched indices of refraction. The details of this new methodology are to be found in:
> A.P. Tran and S.L. Jacques, 2020. Modeling voxel-based Monte Carlo light transport with curved and oblique boundary surfaces. Journal of Biomedical Optics, 25(2), p.025001.

While the information about mcxyz.c are to be found at:
> S. Jacques, T. Li, and S. Prahl, “mcxyz.c, a 3D Monte Carlo simulation of heterogeneous
tissues,” July 2019 version, http://omlc.org/software/mc/mcxyz/index.html .

Similarly to *mcxyz.c*, the extension mcxyzn is written in ANSI standard C, but adapted to handle two forms of parallelism:
> a multi-core CPU implementation using OpenMP that can be run across all platforms (MacOS, LinuxOS, and Windows)
> a GPU implementation that can be used with NVIDIA graphics cards

The 3D Monte Carlo generates an output file of relative fluence rate, F(y,x,z) [W/cm<sup>2</sup> per W delivered] or [1/cm<sup>2</sup>]. The spatial distribution of absorption is obtained by the product of the fluence rate and the absorption coefficient: A(y,x,z) [1/cm<sup>3</sup>] = F(y,x,z) x muav(T(y,x,z), where muav(i) [cm<sup>-1</sup>] is the absorption coefficient of the ith tissue type (the v is for voxel).

