function [figidx] = plot3Dvector(vec,r,xcut,ycut,zcut,figidx,scale)
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

[L,M,N,~] = size(r);

% find indexes for the cuts
[~, xidx] = min(abs(r(:,1,1,1) - xcut));
[~, yidx] = min(abs(r(1,:,1,2) - ycut));
[~, zidx] = min(abs(r(1,1,:,3) - zcut));


% plot figures
figidx = figidx + 1;
figure(figidx);
if (~isempty(scale))
    Scaleplot = scale*ones(M,N);
    imagesc(rot90(Scaleplot));
    hold on;
end
imagesc(rot90(squeeze(vec(xidx,end:-1:1,:))));
set(gca, 'XTick', []);
set(gca, 'YTick', []);
title('Sagittal View of real part');
axis image

figidx = figidx + 1;
figure(figidx);
if (~isempty(scale))
    Scaleplot = scale*ones(L,N);
    imagesc(rot90(Scaleplot));
    hold on;
end
imagesc(rot90(squeeze(vec(:,yidx,:))));
set(gca, 'XTick', []);
set(gca, 'YTick', []);
title('Coronal View of real part');
axis image


figidx = figidx + 1;
figure(figidx);
if (~isempty(scale))
    Scaleplot = scale*ones(L,M);
    imagesc(rot90(Scaleplot));
    hold on;
end
imagesc(rot90(squeeze(vec(:,:,zidx))));
set(gca, 'XTick', []);
set(gca, 'YTick', []);
title('Axial View');
axis image




