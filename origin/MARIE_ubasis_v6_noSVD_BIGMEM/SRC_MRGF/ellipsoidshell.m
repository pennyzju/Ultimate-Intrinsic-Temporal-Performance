function [idx] = ellipsoidshell(r,innR,extR,minZ,maxZ)
%
% _________________________________________________________________________
% _________________________________________________________________________
%
%   Function that generates a region in a cylindric shell
%
% _________________________________________________________________________
% _________________________________________________________________________
%
%

% cylinder = @(r)(abs(r(:,:,:,2) + 1j*((r(:,:,:,1)-0.005)/0.8)) >= innR) & (abs(r(:,:,:,2) + 1j*((r(:,:,:,1)-0.005)/0.8)) <= extR) & (r(:,:,:,3) <= maxZ) & (r(:,:,:,3) >= minZ);
cylinder = @(r)(abs(r(:,:,:,2)-0.005 + 1j*((r(:,:,:,1))/0.85)) >= innR) & (abs(r(:,:,:,2)-0.005 + 1j*((r(:,:,:,1))/0.85)) <= extR) & (r(:,:,:,3) <= maxZ) & (r(:,:,:,3) >= minZ);
pointsI= cylinder(r); % ones for domain elements in the shell
idx = find(pointsI(:)); % get indexes of elements
