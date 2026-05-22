function volRotextended = matrixTranslate(matrix)
   [~, ~, nzmatrix] = size(matrix);
   subplot(2,2,1); imagesc(squeeze(matrix(ceil(end/2),:,:))); axis image off; title('matrix');
   matrixTranslated = imtranslate(matrix,[0,0,nzmatrix],'OutputView','full');
  
   subplot(2,2,2); imagesc(squeeze(matrixTranslated(ceil(end/2),:,:))); axis image off; title('matrixTranslated');
   volRotloose = imrotate3(matrixTranslated,5,[0 1 0],"linear","loose");
   subplot(2,2,3); imagesc(squeeze(volRotloose(ceil(end/2),:,:))); axis image off; title('volRotloose');
   [row,col]  = find(squeeze(volRotloose(ceil(end/2),:,:)));
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
   volRotextended = zeros(nx,ny,nz);
   volRotextended(   xadd/2+1:nxT+xadd/2 , zadd/2+1:nyT+zadd/2 , zadd/2+1:nzT+zadd/2 )      = volRotcrop;
   subplot(2,2,1); imagesc(squeeze(volRotextended(ceil(end/2),:,:))); axis image off; title('volRotextended');

end