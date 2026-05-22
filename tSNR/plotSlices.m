function plotSlices(filepath, title_name, flag)
    % Load the data
    load(filepath);
    
    % Set slice index based on the flag
    switch flag
        case 'front'
            slice_idx = 56;
            get_slice = @(x) squeeze(x(:,slice_idx,:));
        case 'left'
            slice_idx = 48;
            get_slice = @(x) squeeze(x(slice_idx,:,:));
        case 'rot'
            slice_idx = 87;
            get_slice = @(x) squeeze(x(:,:,slice_idx));
        otherwise
            error('Invalid flag. Must be "front", "left", or "rot".');
    end
    
    % Define the output directory
    output_dir = '/data2/jiaxinli/pan/figureplot';
    
    % Create first figure: absolute value plot
    % fig = figure('visible', 'off');  
    % imagesc(abs(get_slice(cv_map(:,:,:))));   
    % axis image off; 
    % title(sprintf('%s %s', title_name, 'alpha'));  % Using input title_name
    % %caxis([0 3e-13]);
    % colormap('jet');
    % colorbar;
    % exportgraphics(fig, fullfile(output_dir, sprintf('%s_alpha.png', title_name)), 'Resolution', 300);
    % close(fig); 

    % Create second figure: log scale plot
    fig = figure('visible', 'off');
    imagesc(log(9*abs(get_slice(mean_map(:,:,:))))); 
    axis image off;
    title(sprintf('%s %s log', title_name, 'alpha'));  % Using input title_name
    caxis([-12 -3]);  
    colormap('jet');
    colorbar;
    % 设置颜色范围
    exportgraphics(fig, fullfile(output_dir, sprintf('%s_log.png', title_name)), 'Resolution', 300);
    close(fig); 
end
