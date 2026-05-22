function [figidx] = plot3Dvector3figs(vec,r,xcut,ycut,zcut,figidx,scale,name,filename)
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

if (~isempty(scale))
    Scaleplot = scale(1,1)*ones(M,N);
    Scaleplot(end,end) = scale(1,2);
    imagesc(rot90(Scaleplot));
    hold on;
end
imagesc(rot90(squeeze(vec(xidx,end:-1:1,:))));
set(gca, 'XTick', []);
set(gca, 'YTick', []);
title(sprintf('(sagittal) %s', name), 'FontSize', 18);
colormap('hot');
colorbar;
axis image
set(gca, 'FontSize', 18);
H = gcf;
saveas(H, sprintf('./FIGS/sagittal_%s',filename), 'epsc');
saveas(H, sprintf('./FIGS/sagittal__%s',filename), 'png');
saveas(H, sprintf('./FIGS/sagittal___%s',filename), 'jpg');
saveas(H, sprintf('./FIGS/sagittal_%s',filename), 'fig');

figidx = figidx + 1;
figure(figidx);
if (~isempty(scale))
    Scaleplot = scale(2,1)*ones(L,N);
    Scaleplot(end,end) = scale(2,2);
    imagesc(rot90(Scaleplot));
    hold on;
end
imagesc(rot90(squeeze(vec(:,yidx,:))));
set(gca, 'XTick', []);
set(gca, 'YTick', []);
title(sprintf('(coronal) %s', name), 'FontSize', 18);
colormap('hot');
colorbar;
axis image
set(gca, 'FontSize', 18);
H = gcf;
saveas(H, sprintf('./FIGS/coronal_%s',filename), 'epsc');
saveas(H, sprintf('./FIGS/coronal__%s',filename), 'png');
saveas(H, sprintf('./FIGS/coronal___%s',filename), 'jpg');
saveas(H, sprintf('./FIGS/coronal_%s',filename), 'fig');

figidx = figidx + 1;
figure(figidx);
if (~isempty(scale))
    Scaleplot = scale(3,1)*ones(L,M);
    Scaleplot(end,end) = scale(3,2);
    imagesc(rot90(Scaleplot));
    hold on;
end
imagesc(rot90(squeeze(vec(:,:,zidx))));
set(gca, 'XTick', []);
set(gca, 'YTick', []);
title(sprintf('(axial) %s', name), 'FontSize', 18);
colormap('hot');
colorbar;
axis image
set(gca, 'FontSize', 18);
H = gcf;
saveas(H, sprintf('./FIGS/axial_%s',filename), 'epsc');
saveas(H, sprintf('./FIGS/axial__%s',filename), 'png');
saveas(H, sprintf('./FIGS/axial___%s',filename), 'jpg');
saveas(H, sprintf('./FIGS/axial_%s',filename), 'fig');



