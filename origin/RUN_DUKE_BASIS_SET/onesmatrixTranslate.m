function onesvolRotextended = onesmatrixTranslate(matrix)
%    matrix = epsilon_r0;
   [~, ~, nzmatrix] = size(matrix);
   subplot(2,2,1); imagesc(squeeze(matrix(ceil(end/2),:,:))); axis image off; title('matrix');
   matrixTranslated = imtranslate(matrix,[0,0,nzmatrix],'OutputView','full');
   matrixTranslated(matrixTranslated<1) = 1;
   subplot(2,2,2); imagesc(squeeze(matrixTranslated(ceil(end/2),:,:))); axis image off; title('matrixTranslated');
   volRotloose = imrotate3(matrixTranslated,5,[0 1 0],"linear","loose");
   volRotloose(volRotloose<1) = 1;
   subplot(2,2,3); imagesc(squeeze(volRotloose(ceil(end/2),:,:))); axis image off; title('volRotloose');
   [row,col]  = find(squeeze(volRotloose(ceil(end/2),:,:))>1);
   [m,index0] = max(row);
   volRotcrop = volRotloose(:,1:m+2,ceil(end/2)-4:ceil(end));
   subplot(2,2,4); imagesc(squeeze(volRotcrop(ceil(end/2),:,:))); axis image off; title('volRotcrop');

   [nxT, nyT, nzT] = size(volRotcrop);
   nx = 95;
   ny = 112;
   nz = 117;
   xadd = nx - nxT;
   %yadd = ny - nyT;
   zadd = nz - nzT;
   onesvolRotextended = ones(nx,ny,nz);
   onesvolRotextended(   xadd/2+1:nxT+xadd/2 , zadd/2+1:nyT+zadd/2 , zadd/2+1:nzT+zadd/2 )      = volRotcrop;
   subplot(2,2,1); imagesc(squeeze(onesvolRotextended(ceil(end/2),:,:))); axis image off; title('onesvolRotextended');

end