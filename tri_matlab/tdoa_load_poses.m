function pose = tdoa_load_poses(file, nb, pos)
% TDOA_LOAD_POSES  Board poses (.R .scale .centre) saved by the calibration in file.
%   Finds the 1 x nb pose struct array (or cell) anywhere in the .mat, whatever its name.
%   pos (optional) = cfg.tri.pos: warns when a saved centre is more than 6 cm away from it
%   (the board was moved: recalibrate).
if nargin < 1 || isempty(file), file = 'tdoa_calibration.mat'; end
if ~isfile(file)
    error('%s not found: calibrate first (tdoa_session, or cfg.freq.calibrate = true)', file);
end
S = load(file);
pose = find_pose(S, nb, 0);
if isempty(pose)
    error('%s: no 1x%d pose (fields R, centre) found. Variables: %s', file, nb, ...
        strjoin(fieldnames(S)', ', '));
end
if nargin >= 3 && ~isempty(pos)
    for b = 1:nb
        p = pos{b}(:)'; if numel(p) == 2, p = [p 0]; end
        cc = pose(b).centre(:)'; if numel(cc) == 2, cc = [cc 0]; end
        if norm(cc(1:2) - p(1:2)) > 0.06
            warning('board %d: calibrated centre is %.0f cm from cfg.tri.pos - recalibrate if it moved', ...
                b, 100 * norm(cc(1:2) - p(1:2)));
        end
    end
end
end

function p = find_pose(S, nb, depth)
p = [];
if depth > 4, return; end
if isstruct(S)
    if isfield(S, 'R') && isfield(S, 'centre') && numel(S) == nb
        p = norm_pose(num2cell(S));
        return
    end
    if numel(S) == 1
        fn = fieldnames(S);
        for k = 1:numel(fn)
            p = find_pose(S.(fn{k}), nb, depth + 1);
            if ~isempty(p), return; end
        end
    end
elseif iscell(S) && numel(S) == nb && all(cellfun(@(x) isstruct(x) && isfield(x, 'R') && isfield(x, 'centre'), S))
    p = norm_pose(S);
end
end

function p = norm_pose(C)
p = struct('R', {}, 'scale', {}, 'centre', {});
for b = 1:numel(C)
    s = 1;
    if isfield(C{b}, 'scale') && ~isempty(C{b}.scale), s = C{b}.scale; end
    cc = C{b}.centre(:)'; if numel(cc) == 2, cc = [cc 0]; end
    p(b) = struct('R', C{b}.R, 'scale', s, 'centre', cc);
end
end
