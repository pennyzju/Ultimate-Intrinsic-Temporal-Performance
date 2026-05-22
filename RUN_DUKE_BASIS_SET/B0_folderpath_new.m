function FolderPath = B0_folderpath_new(B0)
    if B0 == 7
        FolderPath = 'H:\UISNR\20240801_UISNR_output';
    elseif B0 == 3
        %FolderPath = 'F:\20241101_3T_UISNR_output';%front
        FolderPath = 'G:\20251201_3T_UISNR_output';%left,rot
    elseif B0 == 1.5
        FolderPath = 'F:\20251201_1p5T_UISNR_output';
    elseif B0 == 10.5
        FolderPath = '/data2/jiaxinli/pan/20250801_10p5T_UISNR_output';
    elseif B0 == 14
        FolderPath = '/data2/jiaxinli/pan/20250801_14T_UISNR_output';  
    else
        error('B0 must be 7, 3, or 1.5');
    end
    %fprintf('%gT文件夹: %s\n', B0, FolderPath);
end
