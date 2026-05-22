function ubasis_prepare(public_out_dir,out_dir, fileName,flag)
    % function ubasis_prepare(out_dir,epsilon_r,sigma_e,rho,r,f,Dmin,Dmax,regions_dipoles,no_dipole_region)
    load(fullfile(public_out_dir, 'default_material_maps.mat'));
    load(fullfile(public_out_dir, fileName));
    
    if strcmp(flag, 'none')
        idxS = find( mask );
    else
        idxS = find( rotated_mask );
    end
    idxN = find( mask2 );
    idxI = find( mask_dip2 );
    idxS_big = find( mask_big );
    
    % save 
    filename = sprintf('%s/object_def.mat',out_dir);
    % %eval( sprintf('save %s f r dx epsilon_r sigma_e rho idxS idxI idxN mask mask2 dx dy dz pvol',filename) );
    variables_to_save = {'idxS', 'idxI', 'idxN','idxS_big'};
    save(filename, variables_to_save{:});


