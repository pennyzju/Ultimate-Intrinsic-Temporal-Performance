
% addpath(genpath('C:\Users\Guerin\CODE\MARIE_ubasis_v6_noSVD_BIGMEM'));
% addpath(genpath('Z:\ULTIMATE_SNR_SAR\CODE\MARIE_ubasis_v6_noSVD_BIGMEM'));

addpath(genpath('F:\UISNR\rotateback_5\MARIE_ubasis_v6_noSVD_BIGMEM')); % path of the output directory containing all results

close all;

f = 298.0320e6;  % frequency in Hz

% define object/material maps
load DUKE_model_3mm_FLOOD_7T.mat; % data of the mask
xlim = xlim_3mm;   %position of the X pixels in meter
ylim = ylim_3mm;   %position of the Y pixels in meter
zlim = zlim_3mm;   %position of the Z pixels in meter

mask0      = mask_3mm;  % 掩模
epsilon_r0 = perm_3mm;  % 介电常数
sigma_e0   = cond_3mm;  % 电导率
rho0       = dens_3mm;  % 密度

npixadd = 32;  % 2*ceil(Dmax/dx) + 4   *** MUST BE EVEN 必须是偶数***

% extend maps %拓展图像
ind = find(mask0>0);
[indx indy indz] = ind2sub( size(mask0),ind );%下标和线性索引之间进行转换
xmin = min(indx);  xmax = max(indx); 
ymin = min(indy);  ymax = max(indy); 
zmin = min(indz);  zmax = max(indz); 

mask0      = mask0(      xmin:xmax,ymin:ymax,zmin:zmax );
epsilon_r0 = epsilon_r0( xmin:xmax,ymin:ymax,zmin:zmax );
sigma_e0   = sigma_e0(   xmin:xmax,ymin:ymax,zmin:zmax );
rho0       = rho0(       xmin:xmax,ymin:ymax,zmin:zmax );

xlim = xlim(xmin:xmax);
ylim = ylim(ymin:ymax);
zlim = zlim(zmin:zmax);

dx = xlim(2) - xlim(1);
dy = ylim(2) - ylim(1);
dz = zlim(2) - zlim(1);

nx = numel(xlim) + npixadd; %n = numel(A) 返回数组 A 中的元素数目
ny = numel(ylim) + npixadd;
nz = numel(zlim) + npixadd;

lx = dx*nx;
ly = dy*ny;
lz = dz*nz;

x = -lx/2 : dx : lx/2-dx;
y = -ly/2 : dy : ly/2-dy;
z = -lz/2 : dz : lz/2-dz;

[xx, yy, zz] = ndgrid(x,y,z);
r = cat(4,xx,yy,zz);  % 4D vector containing the voxel positions (x = r(:,:,:,1), y = r(:,:,:,2), z = r(:,:,:,3))

% extend material maps to place dipoles 
[nx0, ny0, nz0] = size(mask0);
mask      = zeros(nx,ny,nz);
epsilon_r = ones(nx,ny,nz);
rho       = zeros(nx,ny,nz);
sigma_e   = zeros(nx,ny,nz);

mask(      npixadd/2+1:nx0+npixadd/2 , npixadd/2+1:ny0+npixadd/2 , npixadd/2+1:nz0+npixadd/2 )      = mask0;
epsilon_r( npixadd/2+1:nx0+npixadd/2 , npixadd/2+1:ny0+npixadd/2 , npixadd/2+1:nz0+npixadd/2 )      = epsilon_r0;
rho(       npixadd/2+1:nx0+npixadd/2 , npixadd/2+1:ny0+npixadd/2 , npixadd/2+1:nz0+npixadd/2 )      = rho0;
sigma_e(   npixadd/2+1:nx0+npixadd/2 , npixadd/2+1:ny0+npixadd/2 , npixadd/2+1:nz0+npixadd/2 )      = sigma_e0;


figure;
subplot(2,2,1); imagesc(squeeze(mask(ceil(end/2),:,:))); axis image off; title('mask');
subplot(2,2,2); imagesc(squeeze(epsilon_r(ceil(end/2),:,:))); axis image off; title('eps');
subplot(2,2,3); imagesc(squeeze(rho(ceil(end/2),:,:))); axis image off; title('rho');
subplot(2,2,4); imagesc(squeeze(sigma_e(ceil(end/2),:,:))); axis image off; title('sig');

regions_dipoles = [1 1 1 1 0 1];  % indicate whether to place dipoles in the -x +x -y +y -z +z regions
                                  % beyond the mask

no_dipole_region = [];


%  COMPUTATION OPTIONS
out_dir = './';

Dmin = 3e-2;  % minimum distance between dipoles and object surface 
Dmax = 3.5e-2;  % maximum distance between dipoles and object surface 

Nexc = 1250;  % # of random excitations

tol = 1e-8;  % convergence tolerance for VIE
maxit = 1e6;  % max # of iteration for VIE  

save ubasis_options.mat out_dir Dmin Dmax Nexc tol maxit;  % this file name should not change and should always be in the current directory



% PREPARE
fprintf('PREPARE ...\n'); %输出显示
dipole_prepare(out_dir,epsilon_r0,sigma_e0,rho0,epsilon_r,sigma_e,rho,r,f,Dmin,Dmax,regions_dipoles,no_dipole_region);



% COMPUTE INCIDENT BASIS
fprintf('COMPUTE INCIDENT BASIS ...\n');
ubasis_comp_inc_basis;

% SOLVE FOR INDUCED CURRENTS (J-VIE)
fprintf('SOLVE ...\n');
ubasis_solve;

% COMPUTE SCATTERED FIELDS
fprintf('COMPUTE SCATTERED FIELDS ...\n');
ubasis_comp_fields;












