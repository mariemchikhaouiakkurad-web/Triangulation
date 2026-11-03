function B = tdoa_boards(cfg)
% TDOA_BOARDS  Mic layout, channels and pairs of the 3 boards for the TDOA method.
%   Pairs are taken INSIDE each board only: the mics of one board share one clock,
%   different boards do not (no sync needed).
c = 343;
T = cfg.tri;
loc = {tri_board_layout(cfg.odascfg{1}, cfg.use_ch{1}), ...
       tri_board_layout(cfg.odascfg{2}, cfg.use_ch{2}), ...
       tri_board_layout(cfg.odascfg{3}, T.b_mics)};     % UMA16B: MIC1, MIC16 (read on ch1, ch15)
% the remaining UMA16B mics sit at the board position given in cfg.tri.pos (the cut board
% has no other centre), so their layout is centred on them
loc{3} = loc{3} - mean(loc{3}, 1);
B = struct('name', cfg.names, 'ch', cfg.use_ch, 'Mloc', loc, 'pairs', [], 'maxlag', [], ...
    'weight', num2cell(cfg.tdoa.board_weight), 'centre', T.pos);
for b = 1:numel(B)
    n = size(B(b).Mloc, 1);
    B(b).pairs = nchoosek(1:n, 2);
    d = B(b).Mloc(B(b).pairs(:, 1), :) - B(b).Mloc(B(b).pairs(:, 2), :);
    B(b).maxlag = sqrt(sum(d .^ 2, 2)) / c;
end
end
