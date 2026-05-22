function [idx] = cylindricshell(r,innR,extR,minZ,maxZ)
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

cylinder = @(r)(abs(r(:,:,:,2) + 1j*(r(:,:,:,1))) >= innR) & (abs(r(:,:,:,2) + 1j*(r(:,:,:,1))) <= extR) & (r(:,:,:,3) <= maxZ) & (r(:,:,:,3) >= minZ);
pointsI= cylinder(r); % ones for domain elements in the shell
idx = find(pointsI(:)); % get indexes of elements
