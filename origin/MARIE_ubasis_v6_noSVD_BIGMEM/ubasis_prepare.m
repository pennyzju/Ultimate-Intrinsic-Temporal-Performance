


function ubasis_prepare(out_dir,epsilon_r,sigma_e,rho,r,f,Dmin,Dmax,regions_dipoles,no_dipole_region)

warning('**** In this code, the Z axis of the object MUST correspond to the B0 direction ****');

dx = r(2,1,1,1) - r(1,1,1,1); %resolution in x,y,z
dy = r(1,2,1,2) - r(1,1,1,2);
dz = r(1,1,2,3) - r(1,1,1,3);
pvol = dx * dy * dz;

if dy~=dy || dz~=dx      %dy~=dz
    error('Only isotropic voxels supported.');
end

ind_check = find( epsilon_r(:)<1 | sigma_e(:)<0  | rho(:)<0  );
if ~isempty(ind_check)
    error('Incorrect material maps!');
end

idxS = find( epsilon_r(:)>1 );

% generate dipoles
filename = sprintf('%s/dist_dipolesObject.mat',out_dir);
mask = double( epsilon_r>1 & sigma_e>0 & rho>0 );
figure;
subplot(2,2,1); imagesc(squeeze(mask(ceil(end/2),:,:))); axis image off; title('mask');
subplot(2,2,2); imagesc(squeeze(epsilon_r(ceil(end/2),:,:))); axis image off; title('eps');
subplot(2,2,3); imagesc(squeeze(rho(ceil(end/2),:,:))); axis image off; title('rho');
subplot(2,2,4); imagesc(squeeze(sigma_e(ceil(end/2),:,:))); axis image off; title('sig');

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
[indx indy indz] = ind2sub( size(mask),ind );
xmin = min(indx);  xmax = max(indx); 
ymin = min(indy);  ymax = max(indy); 
zmin = min(indz);  zmax = max(indz); 

[nx ny nz] = size(mask);

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

idxI = find( mask_dip2 );

xd = r(:,:,:,1);
yd = r(:,:,:,2);
zd = r(:,:,:,3);

% position of non-zero voxels
xs = xd(idxS);  
ys = yd(idxS);
zs = zd(idxS);

% position of dipoles
xi = xd(idxI);
yi = yd(idxI);
zi = zd(idxI);

% shave off boundary object voxels
mask2 = ubasis_shave_boundary_voxels(mask);

idxN = find(mask2);

epsN = epsilon_r .* mask2;
sigN = sigma_e .* mask2;
rhoN = rho .* mask2;

mask_dip2 = mask_dip2 + 2*mask + mask2;
ff = figure;
subplot(2,2,1); imagesc(squeeze(mask_dip2(:,:,ceil(end/2)))); axis image off; title('Dipoles, non-zeros and boundary voxels')
subplot(2,2,2); imagesc(squeeze(mask_dip2(:,ceil(end/2),:))); axis image off;
subplot(2,2,3); imagesc(squeeze(mask_dip2(ceil(end/2),:,:))); axis image off; title(sprintf('%d E & B dipole excitations',numel(idxI)))

save mask_dip2.mat mask_dip2;

% save 
filename = sprintf('%s/object_def.mat',out_dir);
eval( sprintf('save %s f r dx epsilon_r sigma_e rho idxS idxI idxN mask mask2 dx dy dz pvol',filename) );




function ubasis_comp_dist_dipoleObject(mask,r,out_dir)


fprintf('Computing object isosurface');
mask_sm = smooth3d(mask,3);
mask_surf = isosurface(r(:,:,:,1),r(:,:,:,2),r(:,:,:,3),mask_sm,0.5);

f=figure;
p=patch(mask_surf,'FaceColor','r'); view(3); axis image; axis off;
lighting gouraud; camlight;
set(p,'FaceColor','red','EdgeColor','red','EdgeLighting','gouraud');
title('Object isosurface'); drawnow;

ind_dip = find(mask==0);
[ind_dip_x ind_dip_y ind_dip_z] = ind2sub(size(mask),ind_dip);

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
end




