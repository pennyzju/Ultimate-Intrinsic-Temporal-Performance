function r = grid3d(x, y, z)

if nargin == 1
    y = x;
    z = x;
end
% [rx, ry, rz] = meshgrid(x, y, z);
% %
% L = length(x);
% M = length(y);
% N = length(z);
% %
% rxx = zeros(L,M,N);
% ryy = zeros(L,M,N);
% %
% r = zeros(L, M, N, 3);
% for index = 1 : N
%     rxx(:,:,index) = rx(:,:,index)';
%     ryy(:,:,index) = ry(:,:,index)';
% end
% r(:,:,:,1) = rxx;
% r(:,:,:,2) = ryy;
% size(rz)
% size(r)
% r(:,:,:,3) = rz;


L = length(x);
M = length(y);
N = length(z);

r = zeros(L,M,N,3);

for ix = 1:L
    xx = x(ix);
    for iy = 1:M
        yy = y(iy);
        for iz = 1:N
            zz = z(iz);
            r(ix,iy,iz,:) = [xx yy zz];
        end
    end
end
% % 
% % % [Xx,Yy,Zz] = meshgrid(x,y,z);
% % % [L,M,N] = size(Xx);
% % % r = zeros(L,M,N,3);
% % % r(:,:,:,1) = Xx;
% % % r(:,:,:,2) = Yy;
% % % r(:,:,:,3) = Zz;