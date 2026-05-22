


function [usnr_map b1_map loss_map cond_number]=compute_usnr_unacc_fast_v2(b1maps,datadir,b1target,mask,lossmatfilepref)
% function [usnr_map b1_map loss_map cond_number]=compute_usnr_unacc_fast_v2(b1maps,datadir,b1target,mask,lossmatfilepref)


[nx ny nz nchannels] = size(b1maps);


% read global SAR matrix
fprintf('\tREADING LOSS MATRIX ...\n');
lossmat=reshape( read_ddata( sprintf('%s/%s.dat',datadir,lossmatfilepref),-1 ),nchannels,2*nchannels );
lossmat=lossmat(:,1:2:end) + 1j*lossmat(:,2:2:end);
q=lossmat;
q_inv=pinv(q.');


ind=find(mask>0);
[indx indy indz]=ind2sub( [nx ny nz],ind );

clear tmp1 tmp2 tmp3;

% delete(gcp);
% parpool( 8 );

% parfor_progress( size(ind,1) );


fprintf('\tSTARTING COMPUTATION ...\n');

% parfor i=1:size(ind,1)

for i=1:size(ind,1)
    
    % parfor_progress();

    fprintf('\t[%d/%d]\n',i,numel(ind));
    
    xpix=indx(i,1);
    ypix=indy(i,1);
    zpix=indz(i,1);
    A=compute_system_matrix(b1maps,xpix,ypix,zpix);
    
    % rf=q_inv*A'*inv( A*q_inv*A' )*b1target;   
    rf=( q_inv*A' )/( A*q_inv*A' )*b1target;  % faster
    
    
    % KKT_matrix = [ q.'  A';
    %                A  0];
               
    tmp4(i,1) = 0;  % cond(KKT_matrix);
    
    b1=A*rf;
    loss_=real( (rf.')*q*conj(rf) );
    
    tmp1(i,1)=abs(b1)/sqrt(loss_);
    tmp2(i,1)=b1;
    tmp3(i,1)=loss_;
    
end

% delete(gcp);
% parfor_progress(0);


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

cond_number = zeros(nx,ny,nz);
cond_number(ind) = tmp4;





        
        
        




