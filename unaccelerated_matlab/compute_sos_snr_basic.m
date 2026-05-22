function [usnr_map, b1_map,loss_map,rf_rate] = compute_sos_snr_basic(B1m, q, mask,idxS)
    % b1maps: 4D array of size [nx, ny, nz, nchannels], containing the B1 maps for each channel
    % q: noise covariance matrix of size [nchannels, nchannels]
    % mask: 3D array of size [nx, ny, nz], indicating the region of interest
    nfields = size(B1m,2);

    nchannels = nfields;
    [nx,ny,nz] = size(mask);
    b1maps = zeros(nx,ny,nz,nchannels);
    tmp = zeros(size(mask));


    for i = 1:nfields
        
        if mod(i,100) == 0
            fprintf('[%d/%d]\n',i,nfields);
        end
        
        tmp( idxS )=B1m(:,i);    
        b1maps(:,:,:,i) = tmp;

    end

    ind = find(mask > 0);  % Find the indices of the mask
    [indx, indy, indz] = ind2sub([nx, ny, nz], ind);

    % Initialize output maps
    b1_map = zeros(nx, ny, nz);
    usnr_map = zeros(nx, ny, nz);

    numInd = size(ind, 1);
    % 预分配内存
    tmp4 = zeros(numInd, 1);
    tmp1 = zeros(numInd, 1);
    tmp2 = zeros(numInd, 1);
    tmp3 = zeros(numInd, 1);

    startTime = tic; 

    % Loop over all pixels in the mask
    for i = 1:numInd
        xpix = indx(i);
        ypix = indy(i);
        zpix = indz(i);

        % Extract the B1 values for all channels at this pixel
        A = compute_system_matrix(b1maps, xpix, ypix, zpix);

        % Compute the root-sum-of-squares (SOS) value
        %b1 = A * A';
        AH = A.';
        b1 = AH'*AH;  

        rf = AH(1:nchannels);      
        loss_=real( (rf')*q*rf);
        %loss_ = real(rf.* q * conj(rf)');
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
end