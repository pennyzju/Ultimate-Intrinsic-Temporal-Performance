function Plot_J(vSolG,vSolG_p,nX,x) 

mSolG = reshape(vSolG, nX, nX, nX, 3);
absJ_VIE = abs3d(mSolG);
mSolG_p = reshape(vSolG_p, nX, nX, nX, 3);
absJ_VIE_p = abs3d(mSolG_p);
%
% figure()
% set(0, 'DefaultTextInterpreter', 'latex');
% subplot(2,3,1)
% imagesc(squeeze(absJ_VIE(:,:,round(nX/2))));
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);
% % title('$z=d$')
% axis image
% subplot(2,3,2)
% imagesc(rot90(squeeze(absJ_VIE(round(nX/2),:,:))));
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);
% % title('$x=d$')
% axis image
% subplot(2,3,3)
% imagesc(x, x, rot90(squeeze(absJ_VIE(:,round(nX/2),:))));
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);
% % title('$y=d$')
% axis image
% subplot(2,3,4)
% imagesc(x, x, squeeze(absJ_VIE_p(:,:,round(nX/2))));
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);
% % title('$z=d$')
% axis image
% subplot(2,3,5)
% imagesc(x, x, rot90(squeeze(absJ_VIE_p(round(nX/2),:,:))));
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);
% % title('$x=d$')
% axis image
% subplot(2,3,6)
% imagesc(x, x, rot90(squeeze(absJ_VIE_p(:,round(nX/2),:))));
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);
% % title('$y=d$')
% axis image
figure()
set(0, 'DefaultTextInterpreter', 'latex');
subplot(2,2,1)
imagesc(squeeze(absJ_VIE(:,:,round(nX/2))));
set(gca, 'XTick', []);
set(gca, 'YTick', []);
% title('$z=d$')
axis image
subplot(2,2,2)
imagesc(rot90(squeeze(absJ_VIE(round(nX/2),:,:))));
set(gca, 'XTick', []);
set(gca, 'YTick', []);
% title('$x=d$')
axis image
subplot(2,2,3)
imagesc(x, x, rot90(squeeze(absJ_VIE(:,round(nX/2),:))));
set(gca, 'XTick', []);
set(gca, 'YTick', []);
% title('$y=d$')
axis image