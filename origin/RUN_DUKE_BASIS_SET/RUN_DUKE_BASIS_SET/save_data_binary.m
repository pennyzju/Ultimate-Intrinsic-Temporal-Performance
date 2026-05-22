

addpath('/autofs/homes/010/guerin/Utilities_matlab/');


load object_def;
load ../basis_ind.mat;
ncomp = size(basis_ind,1);  % # of u-metrics computations


fprintf('DUMPTING GSAR AND LOSS MATRICES -- WITH OUTER VOXEL SHAVED ...\n');
load BASIS_GSAR_N;
load BASIS_LOSS_N;

for i = 1:ncomp

    n=basis_ind(i,1);    
    fprintf('[%d/%d]  N=%d\n',i,ncomp,n);

    LOSS_tw = zeros(n,2*n);    
    LOSS_tw(:,1:2:end) = real( LOSS(1:n,1:n) );
    LOSS_tw(:,2:2:end) = imag( LOSS(1:n,1:n) );    
    
    write_ddata( LOSS_tw,sprintf('LOSS_N_%d.dat',n) );    
    % dlmwrite( sprintf('tot_loss_mat_%d.txt',n),LOSS_tw,'\t' );
    
    GSAR_tw = zeros(n,2*n);
    GSAR_tw(:,1:2:end) = real( GSAR(1:n,1:n) );
    GSAR_tw(:,2:2:end) = imag( GSAR(1:n,1:n) );    
    
    write_ddata(GSAR_tw,sprintf('GSAR_N_%d.dat',n));
    % dlmwrite( sprintf('glob_sar_mat_%d.txt',n),GSAR_tw,'\t' );
    
end

clear GSAR LOSS;



fprintf('DUMPING B1- ...\n');

ind_nzp = find(mask>0);
[indx indy indz] = ind2sub( size(mask),ind_nzp );

xmin = min(indx); xmax = max(indx);  % object' bounding box so we don't store more data than needed
ymin = min(indy); ymax = max(indy);
zmin = min(indz); zmax = max(indz);

nx2 = xmax-xmin+1;
ny2 = ymax-ymin+1;
nz2 = zmax-zmin+1;

write_fdata( epsilon_r( xmin:xmax,ymin:ymax,zmin:zmax ),sprintf('../eps_%d_%d_%d',nx2,ny2,nz2) );
write_fdata( sigma_e( xmin:xmax,ymin:ymax,zmin:zmax ),sprintf('../sig_%d_%d_%d',nx2,ny2,nz2) );
write_fdata( rho( xmin:xmax,ymin:ymax,zmin:zmax ),sprintf('../rho_%d_%d_%d',nx2,ny2,nz2) );

tmp=zeros( size(mask) );

load BASIS_B1m;

nfields = size(B1m,2);
for i = 1:nfields
    
    if mod(i,100) == 0
        fprintf('[%d/%d]\n',i,nfields);
    end
    
    tmp( idxS )=B1m(:,i);    
    
    write_fdata( abs(tmp( xmin:xmax,ymin:ymax,zmin:zmax )),sprintf('B1minus_rowa_coil%d_%d_%d_%d_mag',i,nx2,ny2,nz2) );
    write_fdata( angle(tmp( xmin:xmax,ymin:ymax,zmin:zmax )),sprintf('B1minus_rowa_coil%d_%d_%d_%d_pha',i,nx2,ny2,nz2) );

end

clear B1m;




fprintf('DUMPING B1+ ...\n');

load BASIS_B1p;

tmp=zeros( size(mask) );
nfields = size(B1p,2);
for i = 1:nfields
    
    if mod(i,100) == 0
        fprintf('[%d/%d]\n',i,nfields);
    end
    
    tmp( idxS )=B1p(:,i);    
    
    write_fdata( abs(tmp( xmin:xmax,ymin:ymax,zmin:zmax )),sprintf('B1plus_rowa_coil%d_%d_%d_%d_mag',i,nx2,ny2,nz2) );
    write_fdata( angle(tmp( xmin:xmax,ymin:ymax,zmin:zmax )),sprintf('B1plus_rowa_coil%d_%d_%d_%d_pha',i,nx2,ny2,nz2) );

end

clear B1p;








































