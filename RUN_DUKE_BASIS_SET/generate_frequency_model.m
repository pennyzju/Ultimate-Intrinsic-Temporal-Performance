function generate_frequency_model(infile)
    %GENERATE_FREQUENCY_MODEL
    % Generate material property maps at multiple field strengths
    % from a 7T reference model.
    %
    % CORRECTED VERSION:
    % 1. Protects air/background from scaling.
    % 2. Ensures geometric consistency (mask) across all frequencies.

    % -------------------------------------------------
    % 1. Define supported field strengths
    % -------------------------------------------------
    keys   = {'1p5T', '3T', '7T', '10p5T', '14T'};
    values = {1.5,    3.0,  7.0,  10.5,    14.0};
    
    B0_map = containers.Map(keys, values);
    
    % -------------------------------------------------
    % 2. Load 7T reference model
    % -------------------------------------------------
    if ~exist(infile, 'file')
        error('Input file not found: %s', infile);
    end
    
    S0 = load(infile);
    fprintf('Loaded reference model: %s\n', infile);
    
    % Sanity check
    req_fields = {'cond_3mm','perm_3mm','dens_3mm'};
    for k = 1:numel(req_fields)
        if ~isfield(S0, req_fields{k})
            error('Missing field "%s" in %s', req_fields{k}, infile);
        end
    end
    
    % -------------------------------------------------
    % 3. Create Tissue Mask (CRITICAL STEP)
    % -------------------------------------------------
    % Define what is "Tissue". 
    % Air is typically perm=1.0 and cond=0.0.
    % We select voxels that have significant material properties.
    mask_tissue = (S0.perm_3mm > 1.00) | (S0.cond_3mm > 0);
    
    fprintf('Tissue mask generated. N_voxels = %d\n', nnz(mask_tissue));

    % -------------------------------------------------
    % 4. Loop over field strengths
    % -------------------------------------------------
    B0_names = B0_map.keys;
    
    for i = 1:numel(B0_names)
    
        B0_str = B0_names{i};
        B0     = B0_map(B0_str);
    
        % Output filename
        % Auto-replace '7T' in the filename with the new field strength
        if contains(infile, '7T')
            outfile = strrep(infile, '7T', B0_str);
        else
            % Fallback if filename doesn't contain '7T'
            [~, name, ext] = fileparts(infile);
            outfile = sprintf('%s_%s%s', name, B0_str, ext);
        end
    
        % Skip if exists (Optional: comment out to force overwrite)
        if exist(outfile, 'file')
            fprintf('[SKIP] %s already exists.\n', outfile);
            continue;
        end
    
        % -------------------------------------------------
        % 5. Frequency scaling
        % -------------------------------------------------
        alpha = 0.5;   % conductivity exponent
        beta  = 0.15;  % permittivity exponent
    
        scale_cond = (B0 / 7.0)^alpha;
        scale_perm = (B0 / 7.0)^(-beta);
    
        % Copy structure to initialize
        S = S0;
        
        % -------------------------------------------------
        % 6. Apply Scaling SAFELY
        % -------------------------------------------------
        % A. Initialize with Reference Values
        new_cond = S0.cond_3mm;
        new_perm = S0.perm_3mm;
        
        % B. Scale ONLY the tissue voxels
        % This ensures Air remains exactly Air.
        new_cond(mask_tissue) = double(S0.cond_3mm(mask_tissue)) * scale_cond;
        new_perm(mask_tissue) = double(S0.perm_3mm(mask_tissue)) * scale_perm;
        
        % C. Force Air Properties (Double Safety)
        % Air must have relative permittivity = 1.0 and conductivity = 0.0
        new_cond(~mask_tissue) = 0.0;
        new_perm(~mask_tissue) = 1.0;
        
        % D. Assign back to structure
        S.cond_3mm = new_cond;
        S.perm_3mm = new_perm;
        % S.dens_3mm remains unchanged (density doesn't scale with freq)

        % -------------------------------------------------
        % 7. Save
        % -------------------------------------------------
        save(outfile, '-struct', 'S', '-v7.3');
        fprintf('[DONE] Generated %s (B0 = %.1f T)\n', outfile, B0);
        fprintf('       Scale Factors: Cond=%.3f, Perm=%.3f\n', scale_cond, scale_perm);
    
    end

end