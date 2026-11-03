function [M, N] = tri_geometry(cfg)
% TRI_GEOMETRY  World mic coordinates M{b} (rows = cfg.use_ch{b}) and board normals N{b},
% built exactly as in tri_tdoa_combine: ODAS cfg geometry + pose fitted to measured holes.
nb = numel(cfg.names);
M = cell(1, nb); N = cell(1, nb);
for b = 1:nb
    loc = odas_geometry(cfg.odascfg{b}, cfg.use_ch{b});
    if isfield(cfg, 'geom') && numel(cfg.geom) >= b && ~isempty(cfg.geom{b})
        G = cfg.geom{b};
        if any(isnan(G.xyz(:))), error('%s: fill cfg.geom{%d}.xyz with your measurements (NaN found)', cfg.names{b}, b); end
        if size(G.xyz, 1) >= 3
            Lm = odas_geometry(cfg.odascfg{b}, G.ch);
            [Rb, tb, rms_mm] = kabsch(Lm, G.xyz);
            M{b} = (Rb * loc')' + tb;
            N{b} = Rb(:, 3)';
            fr = [0 0 1];
            if isfield(cfg, 'front') && numel(cfg.front) >= b && ~isempty(cfg.front{b}), fr = cfg.front{b}; end
            if dot(N{b}, fr) < 0, N{b} = -N{b}; end
            fprintf('%-7s pose fitted from %d measured mics, fit error %.1f mm', cfg.names{b}, size(G.xyz, 1), rms_mm);
            if rms_mm > 5, fprintf('  <-- > 5 mm: re-check the measurements'); end
            fprintf('\n');
        else
            M{b} = zeros(numel(cfg.use_ch{b}), 3);
            for i = 1:numel(cfg.use_ch{b})
                r = find(G.ch == cfg.use_ch{b}(i), 1);
                if isempty(r), error('%s: channel %d has no measured position', cfg.names{b}, cfg.use_ch{b}(i)); end
                M{b}(i, :) = G.xyz(r, :);
            end
            N{b} = [0 0 1];
            fprintf('%-7s %d measured mic positions used directly\n', cfg.names{b}, size(M{b}, 1));
        end
    else
        Rb = board_rotation(cfg.facing(b), cfg.roll(b), cfg.mount{b});
        M{b} = (Rb * loc')' + cfg.pos{b};
        N{b} = Rb(:, 3)';
    end
end
end

function loc = odas_geometry(file, ch)
txt = fileread(file);
tok = regexp(txt, 'mu\s*=\s*\(\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*\)', 'tokens');
P = cellfun(@(t) str2double(t), tok, 'UniformOutput', false);
P = vertcat(P{:});
if isempty(P), error('%s: no "mu = (x, y, z)" mic positions found', file); end
mp = regexp(txt, 'map\s*:\s*\(([^)]*)\)', 'tokens', 'once');
if ~isempty(mp)
    map = str2num(mp{1}); %#ok<ST2NM>
else
    map = 1:size(P, 1);
end
loc = zeros(numel(ch), 3);
for i = 1:numel(ch)
    r = find(map == ch(i), 1);
    if isempty(r), error('%s: channel %d not in the mic list', file, ch(i)); end
    loc(i, :) = P(r, :);
end
end

function R = board_rotation(facing, roll, mount)
n = [sind(facing), cosd(facing), 0];
if strcmp(mount, 'flat')
    r = [cosd(facing), -sind(facing), 0]; up = cross([0 0 1], r);
    R = [r', up', [0; 0; 1]];
else
    up = [0 0 1]; r = cross(up, n);
    R = [r', up', n'];
end
Rr = [cosd(roll) -sind(roll) 0; sind(roll) cosd(roll) 0; 0 0 1];
R = R * Rr;
end

function [R, t, rms_mm] = kabsch(L, W)
cl = mean(L, 1); cw = mean(W, 1);
H = (L - cl)' * (W - cw);
[U, ~, V] = svd(H);
D = diag([1 1 sign(det(V * U'))]);
R = V * D * U';
t = cw - (R * cl')';
E = (R * L')' + t - W;
rms_mm = 1000 * sqrt(mean(sum(E .^ 2, 2)));
end
