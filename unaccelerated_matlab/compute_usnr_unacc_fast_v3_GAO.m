

function [usnr_map,b1_map,loss_map,rf_rate,g_factor]=compute_usnr_unacc_fast_v3_GAO(B1m,LOSS,mask,idxS)
% function [usnr_map b1_map loss_map]=compute_usnr_unacc_fast_v2(nx,ny,nz,nchannels,datadir,b1target,mask,lossmatfilepref)


b1target = 1;
% % read B1- map
% b1maps=read_b1maps(datadir,1,nchannels,nx,ny,nz);

[nx,ny,nz] = size(mask);
tmp=zeros( size(mask) );
nfields = size(B1m,2);
nchannels = nfields;
b1maps = zeros(nx,ny,nz,nchannels);
for i = 1:nfields
    
    if mod(i,100) == 0
        fprintf('[%d/%d]\n',i,nfields);
    end
    
    tmp( idxS )=B1m(:,i);    
    b1maps(:,:,:,i) = tmp;
end



% read global SAR matrix

% lossmat=reshape( read_ddata( sprintf('%s/%s.dat',datadir,lossmatfilepref),-1 ),nchannels,2*nchannels );
% lossmat=lossmat(:,1:2:end) + 1j*lossmat(:,2:2:end);
q=LOSS;

ind=find(mask>0);
[indx,indy,indz]=ind2sub( [nx ny nz],ind );

clear tmp1 tmp2 tmp3 tmp4 tmp5;

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
    tmp4(i) =  abs(sum(abs(rf(1:2:length(rf)-1))))./abs(sum(abs(rf(2:2:length(rf)))));
    
    b1=A*rf;
    loss_=real( (rf.')*q*conj(rf) );
    
    tmp1(i,1)=abs(b1)/sqrt(loss_);
    tmp2(i,1)=b1;
    tmp3(i,1)=loss_;
    tmp5(i,1)=sqrt(inv(conj(A)* (q \ conj(A')))*(conj(A)* (q \ conj(A'))));
    % 每1000次输出时间和进度
    if mod(i, 100) == 0 || i == numInd
        elapsedTime = toc(startTime); % 计算已用时间
        fprintf('Processed %d of %d (%.2f%%) - Elapsed time: %.2f seconds\n', ...
                i, numInd, (i / numInd) * 100, elapsedTime);
    end
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

if nargout>3
    rf_rate = zeros(nx,ny,nz);
    rf_rate(ind) = tmp4;
end

if nargout>4
    g_factor = zeros(nx,ny,nz);
    g_factor(ind) = tmp5;
end

% matlabpool close;

% fprintf('\n');




        
        
        




