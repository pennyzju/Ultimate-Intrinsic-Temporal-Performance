function b1maps=read_b1maps(dirpath,nrows,ncoilspr,nx,ny,nz)
% fprintf('Reading B1- maps ...\n');

ind2char='abcdef';
b1maps=zeros( nx,ny,nz,nrows*ncoilspr );

for row=1:nrows

    for chn=1:ncoilspr
        
%         if mod(chn,10)==0
%             fprintf('\tRow #%d out of %d, coil #%d out of %d ...\n',row,nrows,chn,ncoilspr);
%         end

        b1maps(:,:,:,(row-1)*ncoilspr+chn)=reshape( read_fdata(sprintf('%s/B1minus_row%c_coil%d_%d_%d_%d_mag',dirpath,ind2char(row),chn,nx,ny,nz),-1) ...
            .* exp(1j*read_fdata(sprintf('%s/B1minus_row%c_coil%d_%d_%d_%d_pha',dirpath,ind2char(row),chn,nx,ny,nz),-1)),nx,ny,nz );
    end
    
end



