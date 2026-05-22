

function mask2 = ubasis_shave_boundary_voxels(mask)
% function mask2 = ubasis_shave_boundary_voxels(mask)

fprintf('Save off voxels on the boundary of object ...\n');

[nx ny nz] = size(mask);
mask3 = zeros(nx+2,ny+2,nz+2);
mask3(2:end-1,2:end-1,2:end-1) = mask;

ind = find(mask3>0);
[indx indy indz] = ind2sub(size(mask3),ind);
npix = size(ind,1);

mask2 = zeros(size(mask3));
for i=1:npix
    x=indx(i);
    y=indy(i);
    z=indz(i);
    
    if  mask3(x-1,y-1,z-1)>0 && mask3(x-1,y-1,z  )>0 && mask3(x-1,y-1,z+1)>0 && ...
        mask3(x-1,y  ,z-1)>0 && mask3(x-1,y  ,z  )>0 && mask3(x-1,y  ,z+1)>0 && ...
        mask3(x-1,y+1,z-1)>0 && mask3(x-1,y+1,z  )>0 && mask3(x-1,y+1,z+1)>0 && ...        
        mask3(x  ,y-1,z-1)>0 && mask3(x  ,y-1,z  )>0 && mask3(x  ,y-1,z+1)>0 && ...
        mask3(x  ,y  ,z-1)>0 && mask3(x  ,y  ,z  )>0 && mask3(x  ,y  ,z+1)>0 && ...
        mask3(x  ,y+1,z-1)>0 && mask3(x  ,y+1,z  )>0 && mask3(x  ,y+1,z+1)>0 && ...
        mask3(x+1,y-1,z-1)>0 && mask3(x+1,y-1,z  )>0 && mask3(x+1,y-1,z+1)>0 && ...
        mask3(x+1,y  ,z-1)>0 && mask3(x+1,y  ,z  )>0 && mask3(x+1,y  ,z+1)>0 && ...
        mask3(x+1,y+1,z-1)>0 && mask3(x+1,y+1,z  )>0 && mask3(x+1,y+1,z+1)>0     
    mask2(x,y,z) = 1.0;        
    end    
end

mask2 = mask2(2:end-1,2:end-1,2:end-1);
















