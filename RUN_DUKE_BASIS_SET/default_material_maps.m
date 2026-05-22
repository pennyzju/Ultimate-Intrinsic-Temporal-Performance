function default_material_maps(out_dir,B0_str)
    fname = fullfile('E:\UISNR\RUN_DUKE_BASIS_SET', sprintf('DUKE_model_3mm_FLOOD_%s.mat', B0_str));
    DUKE = load(fname);
    
    xlim = DUKE.xlim_3mm;
    ylim = DUKE.ylim_3mm;
    zlim = DUKE.zlim_3mm;

    mask0      = DUKE.mask_3mm;
    epsilon_r0 = DUKE.perm_3mm;
    sigma_e0   = DUKE.cond_3mm;
    rho0       = DUKE.dens_3mm;

    npixadd = 32;  % 2*ceil(Dmax/dx) + 4   *** MUST BE EVEN ***

    % extend maps
    ind = find(mask0>0);
    [indx, indy, indz] = ind2sub( size(mask0),ind );
    xmin = min(indx);  xmax = max(indx); 
    ymin = min(indy);  ymax = max(indy); 
    zmin = min(indz);  zmax = max(indz); 

    mask0      = mask0(      xmin:xmax,ymin:ymax,zmin:zmax );
    epsilon_r0 = epsilon_r0( xmin:xmax,ymin:ymax,zmin:zmax );
    sigma_e0   = sigma_e0(   xmin:xmax,ymin:ymax,zmin:zmax );
    rho0       = rho0(       xmin:xmax,ymin:ymax,zmin:zmax );

    xlim = xlim(xmin:xmax);
    ylim = ylim(ymin:ymax);
    zlim = zlim(zmin:zmax);

    dx = xlim(2) - xlim(1);
    dy = ylim(2) - ylim(1);
    dz = zlim(2) - zlim(1);
    pvol = dx * dy * dz;

    nx = numel(xlim) + npixadd;
    ny = numel(ylim) + npixadd;
    nz = numel(zlim) + npixadd;

    lx = dx*nx;
    ly = dy*ny;
    lz = dz*nz;

    x = -lx/2 : dx : lx/2-dx;
    y = -ly/2 : dy : ly/2-dy;
    z = -lz/2 : dz : lz/2-dz;

    [xx, yy, zz] = ndgrid(x,y,z);
    r = cat(4,xx,yy,zz);  % 4D vector containing the voxel positions (x = r(:,:,:,1), y = r(:,:,:,2), z = r(:,:,:,3))

    % extend material maps to place dipoles 
    [nx0, ny0, nz0] = size(mask0);
    mask      = zeros(nx,ny,nz);
    epsilon_r = ones(nx,ny,nz);
    rho       = zeros(nx,ny,nz);
    sigma_e   = zeros(nx,ny,nz);

    mask(      npixadd/2+1:nx0+npixadd/2 , npixadd/2+1:ny0+npixadd/2 , npixadd/2+1:nz0+npixadd/2 )      = mask0;
    epsilon_r( npixadd/2+1:nx0+npixadd/2 , npixadd/2+1:ny0+npixadd/2 , npixadd/2+1:nz0+npixadd/2 )      = epsilon_r0;
    rho(       npixadd/2+1:nx0+npixadd/2 , npixadd/2+1:ny0+npixadd/2 , npixadd/2+1:nz0+npixadd/2 )      = rho0;
    sigma_e(   npixadd/2+1:nx0+npixadd/2 , npixadd/2+1:ny0+npixadd/2 , npixadd/2+1:nz0+npixadd/2 )      = sigma_e0;


    % 假设 dx, dy, dz 是已定义的浮点数,1e-5定义容差值,避免浮点数比较时由于精度问题引起的误差
    if abs(dy - dx) > 1e-5 || abs(dz - dx) > 1e-5
        error('Only isotropic voxels supported.');
    end

    ind_check = find( epsilon_r(:)<1 | sigma_e(:)<0  | rho(:)<0, 1  );
    if ~isempty(ind_check)
        error('Incorrect material maps!');
    end

    figure;
    subplot(2,2,1); imagesc(squeeze(mask(ceil(end/2),:,:))); axis image off; title('mask');
    subplot(2,2,2); imagesc(squeeze(epsilon_r(ceil(end/2),:,:))); axis image off; title('eps');
    subplot(2,2,3); imagesc(squeeze(rho(ceil(end/2),:,:))); axis image off; title('rho');
    subplot(2,2,4); imagesc(squeeze(sigma_e(ceil(end/2),:,:))); axis image off; title('sig');

    mask = double( epsilon_r>1 & sigma_e>0 & rho>0 );
    mask2 = ubasis_shave_boundary_voxels(mask);

    %%
    load(fullfile(out_dir, 'ubasis_options.mat'))

    % generate dipoles
    filename = sprintf('%s/dist_dipolesObject.mat',out_dir);
    % mask = double( epsilon_r>1 & sigma_e>0 & rho>0 );
    if ~exist(filename,'file')
        ubasis_comp_dist_dipoleObject(mask,r,out_dir);
    else
        load(filename);
        if size(mask_dip) ~= size(mask)
            ubasis_comp_dist_dipoleObject(mask,r,out_dir);
        end
    end
    load(filename);

    % create "no dipole" regions
    ind = find(mask>0);
    [indx, indy, indz] = ind2sub( size(mask),ind );
    xmin = min(indx);  xmax = max(indx); 
    ymin = min(indy);  ymax = max(indy); 
    zmin = min(indz);  zmax = max(indz); 

    [nx, ny, nz] = size(mask);

    if regions_dipoles(1) == 0 && xmin>1  % no dipole in the -x region
        mask_dip( 1:xmin-1,:,: ) = inf;
    end
    if regions_dipoles(2) == 0 && xmax<nx  % +x
        mask_dip( xmax+1:end,:,: ) = inf;
    end
        
    if regions_dipoles(3) == 0 && ymin>1  % no dipole in the -y region
        mask_dip( :,1:ymin-1,: ) = inf;
    end
    if regions_dipoles(4) == 0 && ymax<ny  % +y
        mask_dip( :,ymax+1:end,: ) = inf;
    end

    if regions_dipoles(5) == 0 && zmin>1  % no dipole in the -z region
        mask_dip( :,:,1:zmin-1 ) = inf;
    end
    if regions_dipoles(6) == 0 && zmax<nz  % +z
        mask_dip( :,:,zmax+1:end ) = inf;
    end

    if ~isempty(no_dipole_region)
        mask_dip( find(no_dipole_region>0) ) = inf;
    end

    mask_dip2 = mask_dip>=Dmin &  mask_dip<=Dmax;

    %% 

    save(fullfile(out_dir, 'default_material_maps.mat'),'epsilon_r','sigma_e','rho','r','mask', 'mask2', 'mask_dip2','dx', 'dy', 'dz', 'pvol');


    % mask_dip2 = mask_dip2 + 2*mask + mask2;
    % figure;
    % subplot(2,2,1); imagesc(squeeze(mask_dip2(:,:,ceil(end/2)))); axis image off; title('Dipoles, non-zeros and boundary voxels')
    % subplot(2,2,2); imagesc(squeeze(mask_dip2(:,ceil(end/2),:))); axis image off;
    % subplot(2,2,3); imagesc(squeeze(mask_dip2(ceil(end/2),:,:))); axis image off; title(sprintf('%d E & B dipole excitations',numel(idxI)))

    idxS = find( mask );
    idxN = find( mask2 );
    idxI = find( mask_dip2 );
    %idxS_big = find( mask_big );
    % save 
    filename = sprintf('%s/object_def.mat',out_dir);
    % %eval( sprintf('save %s f r dx epsilon_r sigma_e rho idxS idxI idxN mask mask2 dx dy dz pvol',filename) );
    variables_to_save = {'idxS', 'idxI', 'idxN'};
    save(filename, variables_to_save{:});
