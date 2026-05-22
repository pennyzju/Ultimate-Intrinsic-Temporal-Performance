function iterplot(x,idxSS,Ldim,Mdim,Ndim,ixcut,iycut,izcut,resvec)
%
% _________________________________________________________________________
% _________________________________________________________________________
%
%   Function to plot axial, sagittal and coronal views of the vector
%
%
%
% _________________________________________________________________________
%
%   Computational Prototyping Group, RLE at MIT
% _________________________________________________________________________
% _________________________________________________________________________
%
%


J = zeros(Ldim,Mdim,Ndim,3);
J(idxSS) = x;

colormap('hot');

subplot(2,2,1);
imagesc(rot90(squeeze( abs(J(ixcut,1:end,:,1)) + abs(J(ixcut,1:end,:,2)) + abs(J(ixcut,1:end,:,3)) )));
axis image
title('|J| Sagittal view');

subplot(2,2,2);
imagesc(rot90(squeeze( abs(J(:,iycut,:,1)) + abs(J(:,iycut,:,2)) + abs(J(:,iycut,:,3)) )));
axis image
title('|J| Coronal view');

subplot(2,2,3);
imagesc(rot90(squeeze( abs(J(:,:,izcut,1)) + abs(J(:,:,izcut,2)) + abs(J(:,:,izcut,3)) )));
axis image
title('|J| Axial view');

subplot(2,2,4);
semilogy(1:length(resvec), resvec, 'k-', 'LineWidth', 4);
axis([0 ceil(length(resvec)/100)*100 1e-6 1e1 ]);
grid on;
title('Relative Residue');

snapnow;
% pause(0.001);




