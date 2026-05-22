function default_config_14(out_dir)
f = 588e6;  % frequency in Hz

regions_dipoles = [1 1 1 1 0 1];  % indicate whether to place dipoles in the -x +x -y +y -z +z regions
                                  % beyond the mask

no_dipole_region = [];

Dmin = 3e-2;  % minimum distance between dipoles and object surface 
Dmax = 3.5e-2;  % maximum distance between dipoles and object surface 

Nexc = 1250;  % # of random excitations

tol = 1e-8;  % convergence tolerance for VIE
maxit = 1e6;  % max # of iteration for VIE

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
save(fullfile(out_dir, 'ubasis_options.mat'), 'f', 'regions_dipoles', 'no_dipole_region', 'Dmin', 'Dmax', 'Nexc', 'tol', 'maxit');
