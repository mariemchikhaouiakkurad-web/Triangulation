function tdoa_picture(cfg, stamp)
% TDOA_PICTURE  Picture of a session in the x, y, z reference frame:
%   - axes X, Y, Z from the origin (UMA-8 centre)
%   - every microphone at its calibrated position (tdoa_calibration.mat)
%   - calibration points: true position, measured position, distance true -> measured
%   - live estimates (from <stamp>_tdoa_live.csv)
%   tdoa_picture(cfg, '20261102_093855')   -> live_sessions/<stamp>_frame.png
here = fileparts(mfilename('fullpath'));
dir0 = fullfile(here, 'live_sessions');
S = cfg.tdoa; fs = cfg.fs; blk = round(cfg.block * fs);
B = tdoa_boards(cfg);
L = load(fullfile(here, 'tdoa_calibration.mat'));
for b = 1:3, B(b).pose = L.pose(b); end
cal = S.cal; if size(cal, 2) == 2, cal(:, 3) = S.src_z; end
nCal = size(cal, 1); per = S.move_sec + S.rec_sec;

% measured position at each calibration point (median TDOAs of its recording window)
files = {'uma16a.wav', 'uma8.wav', 'uma16b.wav'};
x = cell(1, 3);
for b = 1:3
    y = audioread(fullfile(dir0, [stamp '_' files{b}]), [1, round(nCal * per * fs)], 'native');
    x{b} = double(y(:, B(b).ch)) / 2 ^ 31;
end
meas = zeros(nCal, 3);
for p = 1:nCal
    tau = cell(1, 3);
    for b = 1:3
        T = [];
        for t0 = (p - 1) * per + S.move_sec + S.settle_sec:0.5:p * per - 0.5
            X = x{b}(round(t0 * fs) + (1:blk), :); X = X - mean(X, 1);
            T(:, end + 1) = tdoa_estimate(X, B(b).pairs, B(b).maxlag, fs, S); %#ok<AGROW>
        end
        tau{b} = median(T, 2);
    end
    meas(p, :) = tdoa_locate(B, tau, cellfun(@(t) ones(size(t)), tau, 'UniformOutput', false), cfg);
end
live = readmatrix(fullfile(dir0, [stamp '_tdoa_live.csv']));
good = live(:, 9) < S.max_rms_us;                    % fit column
lp = live(good, 3:5);

cols = {[0.12 0.47 0.71], [0.93 0.5 0.05], [0.55 0.35 0.75]};
fig = figure('Color', 'w', 'Position', [40 40 1500 820]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Session %s  -  reference frame: origin = UMA-8 centre, X right, Y toward UMA-8, Z up', stamp), ...
    'FontWeight', 'bold', 'FontSize', 14, 'Interpreter', 'none');

% ---------------------------------------------------------------- 3-D frame (big)
a = nexttile(tl, 1, [2 2]); hold(a, 'on'); grid(a, 'on'); axis(a, 'equal'); box(a, 'on');
L1 = 1.0;
quiver3(a, 0, 0, 0, L1, 0, 0, 0, 'Color', [0.85 0.1 0.1], 'LineWidth', 2.5, 'MaxHeadSize', 0.25);
quiver3(a, 0, 0, 0, 0, L1, 0, 0, 'Color', [0.1 0.6 0.1], 'LineWidth', 2.5, 'MaxHeadSize', 0.25);
quiver3(a, 0, 0, 0, 0, 0, L1, 0, 'Color', [0.1 0.25 0.85], 'LineWidth', 2.5, 'MaxHeadSize', 0.25);
text(a, L1 + 0.05, 0, 0, 'X (right)', 'Color', [0.85 0.1 0.1], 'FontWeight', 'bold', 'FontSize', 12);
text(a, 0, L1 + 0.05, 0, 'Y (toward UMA-8)', 'Color', [0.1 0.6 0.1], 'FontWeight', 'bold', 'FontSize', 12);
text(a, 0, 0, L1 + 0.06, 'Z (up)', 'Color', [0.1 0.25 0.85], 'FontWeight', 'bold', 'FontSize', 12);
th = linspace(0, 2 * pi, 120);
for r = [0.5 1]
    plot3(a, r * cos(th), r * sin(th), zeros(size(th)), ':', 'Color', [0.55 0.55 0.55]);
    text(a, r * cosd(-60), r * sind(-60), 0, sprintf('%.1f m', r), 'Color', [0.45 0.45 0.45]);
end
hb = gobjects(1, 3);
for b = 1:3
    W = tdoa_pose(B(b).Mloc, B(b).pose);
    hb(b) = plot3(a, W(:, 1), W(:, 2), W(:, 3), 'o', 'MarkerSize', 5, 'MarkerFaceColor', cols{b}, 'MarkerEdgeColor', 'none');
    c = mean(W, 1);
    text(a, c(1), c(2), c(3) + 0.08, B(b).name, 'Color', cols{b}, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end
hl = plot3(a, lp(:, 1), lp(:, 2), lp(:, 3), '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 8);
for p = 1:nCal
    plot3(a, [cal(p, 1) meas(p, 1)], [cal(p, 2) meas(p, 2)], [cal(p, 3) meas(p, 3)], '--', 'Color', [0.8 0.15 0.12]);
    plot3(a, [0 cal(p, 1)], [0 cal(p, 2)], [0 cal(p, 3)], '-', 'Color', [0.75 0.75 0.75]);
end
ht = plot3(a, cal(:, 1), cal(:, 2), cal(:, 3), 'kx', 'MarkerSize', 13, 'LineWidth', 2.5);
hm = plot3(a, meas(:, 1), meas(:, 2), meas(:, 3), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0.85 0.15 0.12], 'MarkerEdgeColor', 'k');
for p = 1:nCal
    text(a, cal(p, 1), cal(p, 2), cal(p, 3) + 0.07, sprintf('cal %d: %.2f m -> %.2f m', p, norm(cal(p, :)), norm(meas(p, :))), ...
        'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end
legend(a, [hb ht hm hl], {'UMA16A mics', 'UMA-8 mics', 'UMA16B mics (MIC1, MIC16)', 'calibration point (true)', ...
    'calibration point (measured)', 'live estimates (good fit)'}, 'Location', 'northeast');
xlabel(a, 'X (m)'); ylabel(a, 'Y (m)'); zlabel(a, 'Z (m)');
xlim(a, [-1.2 1.2]); ylim(a, [-0.8 1.2]); zlim(a, [-0.2 0.8]);
view(a, -38, 24);
title(a, '3-D reference frame, microphones and calibration points');

% ---------------------------------------------------------------- layout as built (top view, cm)
% Boards from cfg.tri.pos: each UMA-16 mic area (13 cm long) lies along the line from the
% UMA-8 centre to its centre; the UMA-8 board has a 4.5 cm radius.  Your measured distances
% are the GAPS between board edges (14 / 15 cm) and between the board bottoms (18.5 cm).
a = nexttile(tl, 3); hold(a, 'on'); grid(a, 'on'); axis(a, 'equal'); box(a, 'on');
P = 100 * vertcat(cfg.tri.pos{:}); P = P(:, 1:2);
r8 = 4.5; half = 6.5;
th = linspace(0, 2 * pi, 90);
fill(a, r8 * cos(th), r8 * sin(th), cols{2}, 'FaceAlpha', 0.15, 'EdgeColor', cols{2}, 'LineWidth', 1.5);
ends = zeros(2, 2, 2); edge8 = zeros(2, 2);
for k = 1:2
    b = 2 * k - 1;                                    % boards 1 (UMA16A) and 3 (UMA16B)
    u = P(b, :) / norm(P(b, :));
    ends(k, :, 1) = P(b, :) - half * u;               % top end (toward the UMA-8)
    ends(k, :, 2) = P(b, :) + half * u;               % bottom end (XMOS board side)
    edge8(k, :) = r8 * u;
    plot(a, [ends(k, 1, 1) ends(k, 1, 2)], [ends(k, 2, 1) ends(k, 2, 2)], '-', 'Color', cols{b}, 'LineWidth', 7);
end
for b = 1:3                                           % calibrated microphones on top
    W = 100 * tdoa_pose(B(b).Mloc, B(b).pose);
    plot(a, W(:, 1), W(:, 2), 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
end
dimline(a, edge8(1, :), squeeze(ends(1, :, 1)), sprintf('%.1f cm', norm(squeeze(ends(1, :, 1)) - edge8(1, :))), [-3 0]);
dimline(a, edge8(2, :), squeeze(ends(2, :, 1)), sprintf('%.1f cm', norm(squeeze(ends(2, :, 1)) - edge8(2, :))), [3 0]);
dimline(a, squeeze(ends(1, :, 2)), squeeze(ends(2, :, 2)), sprintf('%.1f cm', norm(squeeze(ends(1, :, 2)) - squeeze(ends(2, :, 2)))), [0 -2.5]);
plot(a, [P(1, 1) P(3, 1)], [P(1, 2) P(3, 2)], ':', 'Color', [0.5 0.5 0.5]);
text(a, mean(P([1 3], 1)), mean(P([1 3], 2)) + 1.5, sprintf('centre-centre %.1f cm', norm(P(1, :) - P(3, :))), ...
    'FontSize', 8, 'Color', [0.4 0.4 0.4], 'HorizontalAlignment', 'center');
text(a, r8 + 1, 0, 'UMA-8 (origin)', 'Color', cols{2}, 'FontWeight', 'bold', 'FontSize', 9);
text(a, P(1, 1) - 4, P(1, 2), 'UMA16A', 'Color', cols{1}, 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'FontSize', 9);
text(a, P(3, 1) + 4, P(3, 2), 'UMA16B', 'Color', cols{3}, 'FontWeight', 'bold', 'FontSize', 9);
plot(a, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 1.5);
xlabel(a, 'X (cm)'); ylabel(a, 'Y (cm)');
title(a, 'Layout as built (cm): your measured gaps; o = calibrated mics');

% ---------------------------------------------------------------- distances achieved
a = nexttile(tl, 6); axis(a, 'off');
rows = cell(nCal + 2, 1);
rows{1} = sprintf('%-6s %-20s %-20s %8s %8s', 'point', 'true (x, y, z) m', 'measured (x, y, z) m', 'dist.', 'error');
rows{2} = repmat('-', 1, 68);
for p = 1:nCal
    rows{p + 2} = sprintf('cal %d  (%+.2f %+.2f %+.2f)   (%+.2f %+.2f %+.2f)   %.2f->%.2f %5.0f cm', p, cal(p, :), meas(p, :), ...
        norm(cal(p, :)), norm(meas(p, :)), 100 * norm(meas(p, :) - cal(p, :)));
end
d = vecnorm(lp, 2, 2);
rows{end + 1} = '';
rows{end + 1} = sprintf('live, good fit: %d of %d estimates', size(lp, 1), size(live, 1));
if ~isempty(d)
    rows{end + 1} = sprintf('distance from origin: median %.2f m, max %.2f m', median(d), max(d));
end
text(a, 0, 1, rows, 'FontName', 'Consolas', 'FontSize', 10, 'VerticalAlignment', 'top', 'Interpreter', 'none');
title(a, 'Distances achieved');

png = fullfile(dir0, [stamp '_frame.png']);
exportgraphics(fig, png, 'Resolution', 110);
fprintf('picture saved to %s\n', png);
end

function dimline(a, p, q, label, off)
% dimension line between p and q, label offset by off (cm)
plot(a, [p(1) q(1)], [p(2) q(2)], 'k-', 'LineWidth', 1);
plot(a, [p(1) q(1)], [p(2) q(2)], 'k.', 'MarkerSize', 12);
m = (p + q) / 2 + off;
text(a, m(1), m(2), label, 'FontWeight', 'bold', 'FontSize', 10, 'BackgroundColor', 'w', ...
    'HorizontalAlignment', 'center', 'Margin', 1);
end
