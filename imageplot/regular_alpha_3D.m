% alpha绘图
clear;
clc
%% 7T Pitch
alpha_data = load('H:\UISNR\20240801_UISNR_output\20251022seq\loss_plain\alpha_left_2500.mat');
figure;
imagesc(log(squeeze(alpha_data.alpha(48,:,:))));
subplot(3,6,9)
imagesc(rot90(log(squeeze(alpha_data.alpha(48,14:100,13:104)))));
%imagesc(rot90(log(squeeze(alpha_data.alpha(48,:,:)))));
axis image off; 
colormap('jet');
caxis([-12 2]); 
cb = colorbar;
cb.Ticks = [-12:2:2];     

alpha_data = load('H:\UISNR\20240801_UISNR_output\20251022seq\loss_plain\alpha_wo_left_2500.mat');
subplot(3,6,12)
imagesc(rot90(log(squeeze(alpha_data.alpha(48,14:100,13:104)))));
axis image off; 
colormap('jet');
caxis([-12 2]); 
%colorbar    


%% 7T Roll
alpha_data = load('H:\UISNR\20240801_UISNR_output\20251022seq\loss_plain\alpha_front_2500.mat');

%figure;
subplot(3,6,15)
imagesc(rot90(log(squeeze(alpha_data.alpha(4:90,56,13:104)))));
axis image off; 
colormap('jet');
caxis([-12 2]); 
%colorbar  ; 

alpha_data = load('H:\UISNR\20240801_UISNR_output\20251022seq\loss_plain\alpha_wo_front_2500.mat');
subplot(3,6,18)
imagesc(rot90(log(squeeze(alpha_data.alpha(4:90,56,13:104)))));
axis image off; 
colormap('jet');
caxis([-12 2]); 
%colorbar  ; 


%% 7T Yaw

alpha_data = load('H:\UISNR\20240801_UISNR_output\20251022seq\loss_plain\alpha_rot_2500.mat');
%figure;
subplot(3,6,3)
imagesc(rot90(log(squeeze(alpha_data.alpha(4:90,13:94,87)))));
axis image off; 
colormap('jet');
caxis([-12 2]); 
%colorbar  ;

alpha_data = load('H:\UISNR\20240801_UISNR_output\20251022seq\loss_plain\alpha_wo_rot_2500.mat');
%figure;
subplot(3,6,6)
imagesc(rot90(log(squeeze(alpha_data.alpha(4:90,13:94,87)))));
axis image off; 
colormap('jet');
caxis([-12 2]); 
%colorbar  ;


%% 3T Pitch
alpha_data = load('G:\20251201_3T_UISNR_output\20251022seq\loss_plain\alpha_left_2500.mat');

%figure;
subplot(3,6,8)
imagesc(rot90(log(squeeze(alpha_data.alpha(48,14:100,13:104)))));
axis image off; 
colormap('jet');
caxis([-12 2]); 
%colorbar  ;   
%title('tSNR Slice at X = 48');

alpha_data = load('G:\20251201_3T_UISNR_output\20251022seq\loss_plain\alpha_wo_left_2500.mat');

%figure;
subplot(3,6,11)
imagesc(rot90(log(squeeze(alpha_data.alpha(48,14:100,13:104)))));
axis image off; 
colormap('jet');
caxis([-12 2]); 
%colorbar  ;   
%title('tSNR Slice at X = 48');

%% 3T Roll

alpha_data = load('G:\20251201_3T_UISNR_output\20251022seq\loss_plain\alpha_front_2500.mat');

%figure;
subplot(3,6,14)
imagesc(rot90(log(squeeze(alpha_data.alpha(4:90,56,13:104)))));
axis image off; 
colormap('jet');
caxis([-12 2]); 
%colorbar  ; 

alpha_data = load('G:\20251201_3T_UISNR_output\20251022seq\loss_plain\alpha_wo_front_2500.mat');

%figure;
subplot(3,6,17)
imagesc(rot90(log(squeeze(alpha_data.alpha(4:90,56,13:104)))));
axis image off; 
colormap('jet');
caxis([-12 2]); 
%colorbar  ; 
%% 3T Yaw

alpha_data = load('G:\20251201_3T_UISNR_output\20251022seq\loss_plain\alpha_rot_2500.mat');

%figure;
subplot(3,6,2)
imagesc(rot90(log(squeeze(alpha_data.alpha(4:90,13:94,87)))));
axis image off; 
colormap('jet');
caxis([-12 2]); 
%colorbar  ;


alpha_data = load('G:\20251201_3T_UISNR_output\20251022seq\loss_plain\alpha_wo_rot_2500.mat');

%figure;
subplot(3,6,5)
imagesc(rot90(log(squeeze(alpha_data.alpha(4:90,13:94,87)))));
axis image off; 
colormap('jet');
caxis([-12 2]); 
%colorbar  ;

%% 1.5T Pitch

alpha_data = load('F:\20251201_1p5T_UISNR_output\20251022seq\loss_plain\alpha_left_2500.mat');

%figure;
subplot(3,6,7)
imagesc(rot90(log(squeeze(alpha_data.alpha(48,14:100,13:104)))));
axis image off; 
colormap('jet');
caxis([-12 2]); 
%colorbar  ;   
%title('tSNR Slice at X = 48');

alpha_data = load('F:\20251201_1p5T_UISNR_output\20251022seq\loss_plain\alpha_wo_left_2500.mat');

%figure;
subplot(3,6,10)
imagesc(rot90(log(squeeze(alpha_data.alpha(48,14:100,13:104)))));
axis image off; 
colormap('jet');
caxis([-12 2]); 
%colorbar  ;   
%title('tSNR Slice at X = 48');

%% 1.5T Roll

alpha_data = load('F:\20251201_1p5T_UISNR_output\20251022seq\loss_plain\alpha_front_2500.mat');

%figure;
subplot(3,6,13)
imagesc(rot90(log(squeeze(alpha_data.alpha(4:90,56,13:104)))));
axis image off; 
colormap('jet');
caxis([-12 2]); 
%colorbar  ;   

alpha_data = load('F:\20251201_1p5T_UISNR_output\20251022seq\loss_plain\alpha_wo_front_2500.mat');

%figure;
subplot(3,6,16)
imagesc(rot90(log(squeeze(alpha_data.alpha(4:90,56,13:104)))));
axis image off; 
colormap('jet');
caxis([-12 2]); 
%colorbar  ;   

%% 1.5T Yaw

alpha_data = load('F:\20251201_1p5T_UISNR_output\20251022seq\loss_plain\alpha_rot_2500.mat');

%figure;
subplot(3,6,1)
imagesc(rot90(log(squeeze(alpha_data.alpha(4:90,13:94,87)))));
axis image off; 
colormap('jet');
caxis([-12 2]); 
%colorbar  ;

alpha_data = load('F:\20251201_1p5T_UISNR_output\20251022seq\loss_plain\alpha_wo_rot_2500.mat');

%figure;
subplot(3,6,4)
imagesc(rot90(log(squeeze(alpha_data.alpha(4:90,13:94,87)))));
axis image off; 
colormap('jet');
caxis([-12 2]); 
%colorbar  ;


