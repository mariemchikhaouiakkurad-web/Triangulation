function loc = tri_board_layout(file, ch)
% TRI_BOARD_LAYOUT  Mic coordinates (board frame, m) of channels ch from an ODAS-style cfg.
txt = fileread(file);
tok = regexp(txt, 'mu\s*=\s*\(\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*\)', 'tokens');
P = cellfun(@(t) str2double(t), tok, 'UniformOutput', false);
P = vertcat(P{:});
if isempty(P), error('%s: no "mu = (x, y, z)" mic positions found (run make_geometry_cfgs)', file); end
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
