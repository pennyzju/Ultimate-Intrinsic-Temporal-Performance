function [figidx] = plot1Dcurrent(vec,component,r,figidx)
%
% _________________________________________________________________________
% _________________________________________________________________________
%
%   Function to plot the component value of the currents along each axis
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

[L,M,N,~] = size(r);

xvalues(1:L,1) = r(:,1,1,1);
yvalues(1:M,1) = r(1,:,1,2);
zvalues(1:N,1) = r(1,1,:,3);

xidx = floor(L/2)+1;
yidx = floor(M/2)+1;
zidx = floor(N/2)+1;

% reshape and turn into absolute value
mvec = reshape(vec, L, M, N, 3);

% get values along the axis
Jonx(1:L,:) = mvec(:,yidx,zidx,:); % components on the x axis
Jony(1:M,:) = mvec(xidx,:,zidx,:); % components on the y axis
Jonz(1:N,:) = mvec(xidx,yidx,:,:); % components on the y axis

% plot figures
figidx = figidx + 1;
figure(figidx);
plot(xvalues,abs(Jonx(:,component)));
title('X component of the Electric Current along x axis');


figidx = figidx + 1;
figure(figidx);
plot(yvalues,abs(Jony(:,component)));
title('X component of the Electric Current along y axis');


figidx = figidx + 1;
figure(figidx);
plot(zvalues,abs(Jonz(:,component)));
title('X component of the Electric Current along z axis');


