
% NOT VERY USEFUL...

function [usnr_map b1_map loss_map]=compute_usnr_unacc_fast_v2(nx,ny,nz,nchannels,datadir,b1target,mask,lossmatfilepref)
% function [usnr_map b1_map loss_map]=compute_usnr_unacc_fast_v2(nx,ny,nz,nchannels,datadir,b1target,mask,lossmatfilepref)

% if matlabpool('size')==0
%      matlabpool open 8;  % max number of workers on the icepuffs
% end

% read B1- map
b1maps=read_b1maps(datadir,1,nchannels,nx,ny,nz);

% read global SAR matrix
lossmat=reshape( read_fdata( sprintf('%s/%s.dat',datadir,lossmatfilepref),-1 ),nchannels,2*nchannels );
lossmat=lossmat(:,1:2:end) + 1j*lossmat(:,2:2:end);
q=lossmat;
tic; q_inv=pinv(q.'); toc

% % make sure that it is PSD, i.e. remove small negative eigenvalues due to
% % machine accuracy
% [v d]=eig(lossmat);
% if min(diag(d))<0
%     lossmat=v*abs(d)/v; 
% end

ind=find(mask>0);
[indx indy indz]=ind2sub( [nx ny nz],ind );

% fprintf('\tComputing uSNR ');
% sbuf=[];

clear tmp1 tmp2 tmp3;

% parfor i=1:size(ind,1)
for i=1:size(ind,1)
    
%     for dum=1:size(sbuf,2)
%         fprintf('\b');
%     end
%     sbuf=sprintf('[%d/%d]',i,size(ind,1));
%     fprintf('%s',sbuf);
    
    xpix=indx(i,1);
    ypix=indy(i,1);
    zpix=indz(i,1);
    A=compute_system_matrix(b1maps,xpix,ypix,zpix);
    
    % rf=q_inv*A'*inv( A*q_inv*A' )*b1target;   
    tic;
    rf=( q_inv*A' )/( A*q_inv*A' )*b1target;  % faster
    toc

    tic
    b1=A*rf;
    loss_=real( (rf.')*q*conj(rf) );
    toc
    
    tmp1(i,1)=abs(b1)/sqrt(loss_);
    tmp2(i,1)=b1;
    tmp3(i,1)=loss_;
    
end

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




        
        
        




