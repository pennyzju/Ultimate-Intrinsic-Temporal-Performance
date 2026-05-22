function [idx] = cubicshell(r,xL,yL,zL)
%
% _________________________________________________________________________
% _________________________________________________________________________
%
%   Function that generates a region in a cubic shell
%
% _________________________________________________________________________
% _________________________________________________________________________
%
%

cubeshell = @(r)( r(:,:,:,1) >= max(xL)  );
pointsI= cubeshell(r); % ones for domain elements in the shell

idx = find(pointsI(:)); % get indexes of elements

cubeshell = @(r)( r(:,:,:,1) <= min(xL)  );
pointsI= cubeshell(r); % ones for domain elements in the shell

idx = [idx; find(pointsI(:))]; % get indexes of elements

cubeshell = @(r)( r(:,:,:,2) >= max(yL)  );
pointsI= cubeshell(r); % ones for domain elements in the shell

idx = [idx; find(pointsI(:))]; % get indexes of elements

cubeshell = @(r)( r(:,:,:,2) <= min(yL)  );
pointsI= cubeshell(r); % ones for domain elements in the shell

idx = [idx; find(pointsI(:))]; % get indexes of elements

cubeshell = @(r)( r(:,:,:,3) >= max(zL)  );
pointsI= cubeshell(r); % ones for domain elements in the shell

idx = [idx; find(pointsI(:))]; % get indexes of elements

cubeshell = @(r)( r(:,:,:,3) <= min(zL)  );
pointsI= cubeshell(r); % ones for domain elements in the shell

idx = [idx; find(pointsI(:))]; % get indexes of elements

idx = unique(idx);