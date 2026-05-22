

function [sos_map,b1_map,loss_map]=compute_usnr_sos(B1m,LOSS,mask,idxS)

    [nx,ny,nz] = size(mask);   
    tmp=zeros( size(mask) );
    nfields = size(B1m,2);
    nchannels = nfields;

    
    
    t1 = tic;
    b1maps = zeros(nx,ny,nz,nchannels);
    for i = 1:nfields
        
        if mod(i,100) == 0
            fprintf('[%d/%d]\n',i,nfields);
        end
        
        tmp( idxS )=B1m(:,i);    
        b1maps(:,:,:,i) = tmp;
    end
    elapsedTime = toc(t1)/60;
    fprintf('Elapsed time: %.4f seconds\n', elapsedTime);   
    clear tmp;

    % read global SAR matrix
    % lossmat=reshape( read_ddata( sprintf('%s/%s.dat',datadir,lossmatfilepref),-1 ),nchannels,2*nchannels );
    % lossmat=lossmat(:,1:2:end) + 1j*lossmat(:,2:2:end);
    q=LOSS;
    q = diag(diag(q));
    ind=find(mask>0);
    [indx,indy,indz]=ind2sub( [nx ny nz],ind );
    
    numInd = size(ind, 1);
    
    startTime = tic; % 记录开始时间
    for i=1:numInd
        xpix = indx(i, 1);
        ypix = indy(i, 1);
        zpix = indz(i, 1);
        %b1=A*rf;
        %loss_=real( (rf.')*q*conj(rf) );

        A=reshape( b1maps(xpix,ypix,zpix,:),1,size(b1maps,4));
        %A=A';
        signalmag = A*A';
        rf = A';
        loss_sos = real( (rf.')*q*conj(rf) );
        %noisepower_sos = real(A' * q * A);
        sos_snr = signalmag ./sqrt (noisepower_sos);

        rf = conj( A(1:nchannels) ); 
        loss_sos2=real( (rf.')*q*conj(rf) );
        
        tmp1(i,1)=sos_snr;
        tmp2(i,1)=signalmag ;
        tmp3(i,1)=noisepower_sos;
        % 每1000次输出时间和进度
        if mod(i, 1000) == 0 || i == numInd
            elapsedTime = toc(startTime)/60; % 计算已用时间
            fprintf('Processed %d of %d (%.2f%%) - Elapsed time: %.2f minutes\n', ...
                    i, numInd, (i / numInd) * 1000, elapsedTime);
        end
    end
       
    sos_map=zeros(nx,ny,nz);
    sos_map(ind)=tmp1;
    
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
    
    
    
    
    
            
            
            
    
    
    
    
    