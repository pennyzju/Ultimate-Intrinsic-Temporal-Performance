


function [usnr_map, b1_map, loss_map]=compute_usnr_unacc_fast_v3(nx,ny,nz,nchannels,datadir,b1target,mask,lossmatfilepref)
% function [usnr_map b1_map loss_map]=compute_usnr_unacc_fast_v2(nx,ny,nz,nchannels,datadir,b1target,mask,lossmatfilepref)
% compute_usnr_unacc_fast_v3(63,80,85,2500,"/data/jiaxinli/projects/coil/results/20240601",1e-3,mask,"LOSS_N_")

if isempty(gcp('nocreate'))
    parpool('local', 8);  % max number of workers on the icepuffs
end


% read B1- map
b1maps=read_b1maps(datadir,1,nchannels,nx,ny,nz);


% read global SAR matrix
lossmat=reshape( read_ddata( sprintf('%s/%s.dat',datadir,lossmatfilepref),-1 ),nchannels,2*nchannels );
lossmat=lossmat(:,1:2:end) + 1j*lossmat(:,2:2:end);
q=lossmat;

ind_nzp = find(mask>0);
[indx, indy, indz] = ind2sub( size(mask),ind_nzp );

xmin = min(indx); xmax = max(indx);  % object' bounding box so we don't store more data than needed
ymin = min(indy); ymax = max(indy);
zmin = min(indz); zmax = max(indz);

nx2 = xmax-xmin+1;
ny2 = ymax-ymin+1;
nz2 = zmax-zmin+1;
mask2 = mask(xmin:xmax,ymin:ymax,zmin:zmax);
ind = find(mask2>0);
[indx, indy, indz] = ind2sub( size(mask2),ind );

clear tmp1 tmp2 tmp3;

parfor_progress(size(ind,1));  % initialize

parfor i=1:size(ind,1)
    
    parfor_progress;  % print progress
    
    xpix = indx(i,1);
    ypix = indy(i,1);
    zpix = indz(i,1);
    A = compute_system_matrix(b1maps,xpix,ypix,zpix);

    A2 = [q A.' ; conj(A) 0];
    b2 = [zeros(nchannels,1) ; b1target ];
    
    % sol = A2 \ b2;
    sol = pinv(A2) * b2;
    
    rf = conj( sol(1:nchannels) );
    
    b1=A*rf;
    loss_=real( (rf.')*q*conj(rf) );
    
    tmp1(i,1)=abs(b1)/sqrt(loss_);
    tmp2(i,1)=b1;
    tmp3(i,1)=loss_;
    
end

parfor_progress(0);  % clean up

usnr_map=zeros(nx,ny,nz);
usnr_map(ind)=tmp1;

if nargout>1
    b1_map=zeros(nx,ny,nz);
    b1_map(ind)=tmp2;
end

if nargout>2
    loss_map=zeros(nx,ny,nz);
    loss_map(ind)=tmp3;
end

% matlabpool close;

% fprintf('\n');




        
        
        




