function [figidx] = plotEps(e_r,r,xcut,ycut,zcut,figidx)
%
% _________________________________________________________________________
% _________________________________________________________________________
%
%   Function to plot axial, sagittal and coronal views of the epsilon
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

if (nargin < 6) || isempty(figidx)
    figidx = 0;
end


% find indexes for the cuts
[~, xidx] = min(abs(r(:,1,1,1) - xcut));
[~, yidx] = min(abs(r(1,:,1,2) - ycut));
[~, zidx] = min(abs(r(1,1,:,3) - zcut));

realeps = real(e_r);
imageps = imag(e_r);

% plot figures
figidx = figidx + 1;
figure(figidx);
imagesc(rot90(squeeze(realeps(xidx,end:-1:1,:))));
set(gca, 'XTick', []);
set(gca, 'YTick', []);
title('Sagittal View of real part');
axis image

figidx = figidx + 1;
figure(figidx);
imagesc(rot90(squeeze(imageps(xidx,end:-1:1,:))));
set(gca, 'XTick', []);
set(gca, 'YTick', []);
title('Sagittal View of imaginary part');
axis image


figidx = figidx + 1;
figure(figidx);
imagesc(rot90(squeeze(realeps(:,yidx,:))));
set(gca, 'XTick', []);
set(gca, 'YTick', []);
title('Coronal View of real part');
axis image

figidx = figidx + 1;
figure(figidx);
imagesc(rot90(squeeze(imageps(:,yidx,:))));
set(gca, 'XTick', []);
set(gca, 'YTick', []);
title('Coronal View of imaginary part');
axis image


figidx = figidx + 1;
figure(figidx);
imagesc(rot90(squeeze(realeps(:,:,zidx))));
set(gca, 'XTick', []);
set(gca, 'YTick', []);
title('Axial View');
axis image

figidx = figidx + 1;
figure(figidx);
imagesc(rot90(squeeze(imageps(:,:,zidx))));
set(gca, 'XTick', []);
set(gca, 'YTick', []);
title('Axial View');
axis image

