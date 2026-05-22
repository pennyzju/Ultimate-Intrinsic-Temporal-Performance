function result = analyze_usnr_convergence(data_dir, voxels, flag, angle)
%ANALYZE_USNR_CONVERGENCE  Paper-level convergence analysis for USNR
%
% result = analyze_usnr_convergence(data_dir, voxels, flag, angle)
%
% result fields:
%   .suffix
%   .usnr        [Nfile x Nvoxel]
%   .rel_change  [Nfile-1 x Nvoxel]

%% File search
pattern = sprintf('USNR_%s_%d_*.mat', flag, angle);
files = dir(fullfile(data_dir, pattern));
assert(~isempty(files), 'No USNR files found.');

%% Extract suffix
Nf = numel(files);
suffix = zeros(Nf,1);
for i = 1:Nf
    tok = regexp(files(i).name, '_(\d+)\.mat$', 'tokens');
    suffix(i) = str2double(tok{1}{1});
end
[suffix, idx] = sort(suffix);
files = files(idx);

Nv = size(voxels,1);
usnr = zeros(Nf, Nv);

%% Load and sample
for i = 1:Nf
    S = load(fullfile(data_dir, files(i).name));
    fn = fieldnames(S);

    % auto-detect 3D matrix
    for k = 1:numel(fn)
        tmp = S.(fn{k});
        if isnumeric(tmp) && ndims(tmp)==3
            U = tmp;
            break;
        end
    end

    for v = 1:Nv
        usnr(i,v) = U(voxels(v,1), voxels(v,2), voxels(v,3));
    end
end

%% Relative change
rel_change = abs(diff(usnr)) ./ abs(usnr(1:end-1,:));

%% Pack result
result.suffix     = suffix;
result.usnr       = usnr;
result.rel_change = rel_change;
result.voxels     = voxels;
result.flag       = flag;
result.angle      = angle;

end
