clearvars -except B1m LOSS
usnr_mask = zeros(95,112,117);
usnr_mask(48,:,:) = 1;
load('H:\20240801_UISNR_output\left_0\left_0_material_maps.mat','rotated_mask');
load('H:\20240801_UISNR_output\left_0\object_def.mat','idxS');
usnr_mask_rot = usnr_mask .* rotated_mask;

[usnr_map_original, b1_map_original, loss_map_original, rf_rate] = compute_sos_snr_basic(B1m, LOSS, usnr_mask_rot, idxS);

coilloss = diag(diag(LOSS));
totLOSS = LOSS + coilloss;
[usnr_map_coilloss, b1_map_coilloss, loss_map_coilloss, rf_rate] = compute_sos_snr_basic(B1m, totLOSS, usnr_mask_rot, idxS);

figure;
clim=[0 1e15];
subplot(1,3,1);
imagesc(squeeze(loss_map_coilloss(48,:,:)));
title('coilloss');
clim;
subplot(1,3,2);
imagesc(squeeze(loss_map_original(48,:,:)));
title('original')
clim;
colorbar
test = loss_map_coilloss./loss_map_original;
figure;
imagesc(squeeze(test(48,:,:)));
colorbar