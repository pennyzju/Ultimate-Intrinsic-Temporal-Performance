clear;

%load('/data2/jiaxinli/pan/20240801_UISNR_output/front_-1/front_-1_material_maps.mat')
%load('/data2/jiaxinli/pan/20240801_UISNR_output/public/default_material_maps.mat')
load('H:\20240801_UISNR_output\left_0\default_material_maps.mat');
load('H:\20240801_UISNR_output\left_0\left_0_material_maps.mat')
filename = 'head_surface' ;
% 构建网格
    [L,M,N] = size(mask);
    [x,y,z] = meshgrid(1:M, 1:L, 1:N);

    % 平滑处理（可选）
    mask_smooth = smooth3(mask);

    % 绘制表面
    figure;
    p = patch(isosurface(x,y,z,mask_smooth,0.5));
    set(p,'FaceColor','red','EdgeColor','none');  

    % 视觉效果
    daspect([1 1 1]);
    %view(3); 
    view([1 1 1])
    camlight; 
    lighting gouraud;
    axis off;
%view([90 0]);   % 把相机转到反方向

    % 保存为 fig 和 png
    savefig([filename '.fig']);  
    saveas(gcf, [filename '.png']);  

idxI = find( mask_dip2 == 1);
N_dipole = numel(idxI);
[L,M,N] = size(mask_dip2);
nD = L*M*N;
idxI3 = [idxI; nD+idxI; 2*nD+idxI]; % the vector of input positions in 3D grid

xd = floor(r(:,:,:,1)./0.002)+40;
yd = floor(r(:,:,:,2)./0.002)+45;
zd = double(int8(r(:,:,:,3)./0.002)+50);

ys = yd(idxI);
zs = zd(idxI);

xs = xd(idxI);  

uiopen('./head_surface.fig',1)
hold on;

[x,y,z]=meshgrid(1:77,1:96,1:87);
kk=(permute(smooth3(mask_dip2==1),[3,1,2]));
p = patch(isosurface(x,y,z,double(kk)))
set(p,'FaceColor',[0.94,0.94,0.94],'LineStyle','none')
set(gcf,'Position',[30 400 800 600]);
set(gcf,'Color','w');

set(gcf,'Position',[30 400 800 600]);
set(gcf,'Color','w');
saveas(gcf,['./monkey_model_currrent_surf.tiff']);