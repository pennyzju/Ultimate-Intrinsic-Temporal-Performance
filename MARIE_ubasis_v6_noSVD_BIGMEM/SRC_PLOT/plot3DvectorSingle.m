function [figidx] = plot3DvectorSingle(vec,r,xcut,ycut,zcut,figidx,scale,name)
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

if (nargin < 6) || isempty(figidx)
    figidx = 0;
end
if (nargin < 7) || isempty(scale)
    scale = [];
end
if (nargin < 8) || isempty(name)
    name = ' ';
end


[L,M,N,~] = size(r);

% find indexes for the cuts
[~, xidx] = min(abs(r(:,1,1,1) - xcut));
[~, yidx] = min(abs(r(1,:,1,2) - ycut));
[~, zidx] = min(abs(r(1,1,:,3) - zcut));


% plot figures
figidx = figidx + 1;
figure(figidx);

subplot(1,3,1);
if (~isempty(scale))
    Scaleplot = scale*ones(M,N);
    imagesc(rot90(Scaleplot));
    hold on;
end
imagesc(rot90(squeeze(vec(xidx,end:-1:1,:))));
set(gca, 'XTick', []);
set(gca, 'YTick', []);
title(sprintf('SAGGITAL %s', name));
colorbar
axis image

subplot(1,3,2);
if (~isempty(scale))
    Scaleplot = scale*ones(L,N);
    imagesc(rot90(Scaleplot));
    hold on;
end
imagesc(rot90(squeeze(vec(:,yidx,:))));
set(gca, 'XTick', []);
set(gca, 'YTick', []);
title(sprintf('CORONAL %s', name));
colorbar
axis image


subplot(1,3,3);
if (~isempty(scale))
    Scaleplot = scale*ones(L,M);
    imagesc(rot90(Scaleplot));
    hold on;
end
imagesc(rot90(squeeze(vec(:,:,zidx))));
set(gca, 'XTick', []);
set(gca, 'YTick', []);
title(sprintf('AXIAL %s', name));
colorbar;
axis image




