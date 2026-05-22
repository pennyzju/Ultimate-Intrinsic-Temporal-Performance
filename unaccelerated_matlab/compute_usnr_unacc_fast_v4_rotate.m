

function [usnr_map,b1_map,loss_map,rf_rate]=compute_usnr_unacc_fast_v4_rotate(B1m,LOSS,mask,idxS,angle_deg,point1,point2)

    b1target = 1e-3;
    % % read B1- map
    % b1maps=read_b1maps(datadir,1,nchannels,nx,ny,nz);

    [nx,ny,nz] = size(mask);

    % 创建一个与 A 同尺寸的零矩阵mask
    mask_zeroed = zeros(nx,ny,nz);
    

    % 只保留中间层，将其他层数据置零
    mask_zeroed(ceil(nx / 2),:,:) = mask(ceil(nx / 2),:,:);
    mask_zeroed(:,ceil(ny / 2),:) = mask(:,ceil(ny / 2),:);
    mask_zeroed(:,:,ceil(nz / 2)) = mask(:,:,ceil(nz / 2));

    mask = mask_zeroed;
    mask = rotate_matrix(mask, point1, point2, angle_deg);

    
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
    
    
    
    clear tmp;
    % read global SAR matrix
    
    % lossmat=reshape( read_ddata( sprintf('%s/%s.dat',datadir,lossmatfilepref),-1 ),nchannels,2*nchannels );
    % lossmat=lossmat(:,1:2:end) + 1j*lossmat(:,2:2:end);
    q=LOSS;
    
    ind=find(mask>0);
    [indx,indy,indz]=ind2sub( [nx ny nz],ind );
    
    numInd = size(ind, 1);
    % 预分配内存，使用 gpuArray 来将数组转移到 GPU
    tmp4 = gpuArray.zeros(numInd, 1);
    tmp1 = gpuArray.zeros(numInd, 1);
    tmp2 = gpuArray.zeros(numInd, 1);
    tmp3 = gpuArray.zeros(numInd, 1);
    
    % 将常量或不变的变量转移到 GPU
    q_gpu = gpuArray(q);
    b1target_gpu = gpuArray(b1target);
    indx_gpu = gpuArray(indx);
    indy_gpu = gpuArray(indy);
    indz_gpu = gpuArray(indz);
    
    
    startTime = tic; % 记录开始时间
    for i=1:numInd
        xpix = indx_gpu(i, 1);
        ypix = indy_gpu(i, 1);
        zpix = indz_gpu(i, 1);
        
        % 将 compute_system_matrix 的输出转移到 GPU
        A = gpuArray(compute_system_matrix(b1maps, xpix, ypix, zpix));
         % 构建 A2 和 b2，利用 GPU 加速
        A2 = [q_gpu, A.'; conj(A), 0];
        b2 = [zeros(nchannels, 1, 'gpuArray'); b1target_gpu];
    
        % sol = A2 \ b2;
        sol = pinv(A2) * b2;
        
        rf = conj( sol(1:nchannels) );  
        tmp4(i) =  abs(sum(abs(rf(1:2:length(rf)-1))))./abs(sum(abs(rf(2:2:length(rf)))));
        
        b1=A*rf;
        loss_=real( (rf.')*q*conj(rf) );
        
        tmp1(i,1)=abs(b1)/sqrt(loss_);
        tmp2(i,1)=b1;
        tmp3(i,1)=loss_;
        % 每1000次输出时间和进度
        if mod(i, 100) == 0 || i == numInd
            elapsedTime = toc(startTime); % 计算已用时间
            fprintf('Processed %d of %d (%.2f%%) - Elapsed time: %.2f seconds\n', ...
                    i, numInd, (i / numInd) * 100, elapsedTime);
        end
    end
    % 将结果从 GPU 传回 CPU
    tmp4 = gather(tmp4);
    tmp1 = gather(tmp1);
    tmp2 = gather(tmp2);
    tmp3 = gather(tmp3);
    
    
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
    
    
    % matlabpool close;
    
    % fprintf('\n');
    
    
    
    
            
            
            
    
    
    
    
    