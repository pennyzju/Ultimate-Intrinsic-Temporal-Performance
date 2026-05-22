function ubasis_comp_dist_dipoleObject(mask,r,out_dir)


fprintf('Computing object isosurface');
mask_sm = smooth3d(mask,3);
mask_surf = isosurface(r(:,:,:,1),r(:,:,:,2),r(:,:,:,3),mask_sm,0.5);

figure;
p=patch(mask_surf,'FaceColor','r'); view(3); axis image; axis off;
lighting gouraud; camlight;
set(p,'FaceColor','red','EdgeColor','red','EdgeLighting','gouraud');
title('Object isosurface'); drawnow;

ind_dip = find(mask==0);
[ind_dip_x, ind_dip_y, ind_dip_z] = ind2sub(size(mask),ind_dip);

nv = size(mask_surf.vertices,1);  % # of vertices

mask_dip = zeros(size(mask));
xs = squeeze( r(:,1,1,1) );
ys = squeeze( r(1,:,1,2) );
zs = squeeze( r(1,1,:,3) );

fprintf('Computing minimum distance between dipole voxels and object surface ');
sbuf = [];
for ii = 1:size(ind_dip,1)
    
    if mod(ii,1000)==0
        for j=1:size(sbuf,2)
            fprintf('\b');
        end
        sbuf = sprintf('[%d/%d]\n',ii,size(ind_dip,1));
        fprintf(sbuf);
    end
    
    % dipole location
    x = xs( ind_dip_x(ii) );
    y = ys( ind_dip_y(ii) );
    z = zs( ind_dip_z(ii) );
    
    % find min distance between voxel and object surface
    dist = mask_surf.vertices - repmat([x y z],[nv 1]);
    dist = sqrt( sum( abs(dist).^2, 2) );
    mask_dip(ind_dip(ii)) = min(dist(:));
end
fprintf('\n');

eval( sprintf('save %s/dist_dipolesObject.mat mask_dip mask_surf mask_sm',out_dir) )
