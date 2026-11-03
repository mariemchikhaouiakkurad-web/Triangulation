function tdoa_evaluate(cfg, stamp, win, ~)
% TDOA_EVALUATE  Calibrate + localise a recorded tdoa_session and measure the error at the
% test points (known positions).  Can be rerun on any saved session:
%   tdoa_evaluate(cfg, '20261102_101500')       % windows rebuilt from cfg.tdoa timings
%   win  (n_cal + n_test) x 2 recording windows (s)
%
% For every test point: true / measured x, y, z, distance from the origin (UMA-8 centre),
% distance error, direction error (deg, seen from the centre of the boards).
% Writes <stamp>_tdoa_eval.csv and the picture <stamp>_tests.png in live_sessions.
S = cfg.tdoa;
cal = S.cal; test = S.test;
if size(cal, 2) == 2, cal(:, 3) = S.src_z; end
if size(test, 2) == 2, test(:, 3) = S.src_z; end
steps = [cal; test];
nCal = size(cal, 1); nTest = size(test, 1);
if nargin < 3 || isempty(win)
    per = S.move_sec + S.rec_sec;
    k = (1:size(steps, 1))';
    win = [(k - 1) * per + S.move_sec + S.settle_sec, k * per];
end
outdir = fullfile(fileparts(mfilename('fullpath')), 'live_sessions');
res = tdoa_offline(cfg, stamp, win(1:nCal, :));
o = mean(vertcat(cfg.tri.pos{:}), 1); if numel(o) == 2, o(3) = 0; end

R = nan(nTest, 8);                                    % meas xyz, dist true, dist meas, dir err, spread, n
fprintf('\n=== accuracy of the own TDOA method at the test points ===\n');
fprintf('%-7s %-21s %-21s %8s %8s %9s %9s %16s\n', 'point', 'true (x, y, z) m', 'measured (x, y, z) m', ...
    'dist', 'meas.', 'dist.err', 'dir.err', 'dx/dy/dz cm');
for s = 1:nTest
    i = nCal + s;
    in = res.t >= win(i, 1) & res.t <= win(i, 2) & all(isfinite(res.pos), 2)';
    if ~any(in), fprintf('test %-2d no loud block (phone not heard)\n', s); continue; end
    m = median(res.pos(in, :), 1);
    vt = test(s, :) - o; vm = m - o;
    direrr = acosd(max(-1, min(1, dot(vt, vm) / (norm(vt) * norm(vm)))));
    sp = median(vecnorm(res.pos(in, :) - m, 2, 2));
    R(s, :) = [m, norm(test(s, :)), norm(m), direrr, sp, sum(in)];
    d = m - test(s, :);
    fprintf('test %-2d (%+.2f %+.2f %+.2f)  (%+.2f %+.2f %+.2f)  %6.2f m %6.2f m %+7.0f cm %7.1f deg  %+4.0f/%+4.0f/%+4.0f\n', ...
        s, test(s, :), m, R(s, 4), R(s, 5), 100 * (R(s, 5) - R(s, 4)), direrr, 100 * d);
end
ok = all(isfinite(R(:, 1:3)), 2);
if any(ok)
    e3 = vecnorm(R(ok, 1:3) - test(ok, :), 2, 2);
    fprintf(['\nsummary over %d test points: median 3-D error %.0f cm, median direction error %.1f deg, ' ...
             'median distance error %+.0f cm\n'], sum(ok), 100 * median(e3), median(R(ok, 6)), 100 * median(R(ok, 5) - R(ok, 4)));
end
ev = fopen(fullfile(outdir, [stamp '_tdoa_eval.csv']), 'w');
fprintf(ev, 'point,true_x,true_y,true_z,meas_x,meas_y,meas_z,true_dist_m,meas_dist_m,dist_err_m,dir_err_deg,spread_m,blocks\n');
for s = 1:nTest
    fprintf(ev, 'test %d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.2f,%.3f,%d\n', s, test(s, :), R(s, 1:3), ...
        R(s, 4), R(s, 5), R(s, 5) - R(s, 4), R(s, 6), R(s, 7), R(s, 8));
end
fclose(ev);

% ---------------------------------------------------------------- picture
fig = figure('Color', 'w', 'Position', [60 60 1400 820], 'Name', ['tests ' stamp]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Test session %s  -  x = true position, o = measured (origin = UMA-8 centre)', stamp), ...
    'FontWeight', 'bold', 'FontSize', 14, 'Interpreter', 'none');
P = vertcat(cfg.tri.pos{:}); if size(P, 2) == 2, P(:, 3) = 0; end
cols = lines(nTest);
% top view
a = nexttile(tl); hold(a, 'on'); grid(a, 'on'); axis(a, 'equal'); box(a, 'on');
plot(a, P(:, 1), P(:, 2), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', [0.5 0.5 0.5]);
plot(a, cal(:, 1), cal(:, 2), '+', 'Color', [0.6 0.6 0.6], 'MarkerSize', 9, 'LineWidth', 1.5);
for s = 1:nTest
    plot(a, test(s, 1), test(s, 2), 'x', 'Color', cols(s, :), 'MarkerSize', 12, 'LineWidth', 2.5);
    text(a, test(s, 1), test(s, 2), sprintf('  %d', s), 'Color', cols(s, :), 'FontWeight', 'bold');
    if ok(s)
        plot(a, [test(s, 1) R(s, 1)], [test(s, 2) R(s, 2)], '--', 'Color', cols(s, :));
        plot(a, R(s, 1), R(s, 2), 'o', 'MarkerSize', 9, 'MarkerFaceColor', cols(s, :), 'MarkerEdgeColor', 'k');
    end
end
xlabel(a, 'X (m)  right'); ylabel(a, 'Y (m)  toward UMA-8'); title(a, 'Top view (+ = calibration points)');
% side view: horizontal distance from the origin vs height
a = nexttile(tl); hold(a, 'on'); grid(a, 'on'); axis(a, 'equal'); box(a, 'on');
yline(a, 0, 'k-');
for s = 1:nTest
    plot(a, hypot(test(s, 1), test(s, 2)), test(s, 3), 'x', 'Color', cols(s, :), 'MarkerSize', 12, 'LineWidth', 2.5);
    if ok(s)
        hm = hypot(R(s, 1), R(s, 2));
        plot(a, [hypot(test(s, 1), test(s, 2)) hm], [test(s, 3) R(s, 3)], '--', 'Color', cols(s, :));
        plot(a, hm, R(s, 3), 'o', 'MarkerSize', 9, 'MarkerFaceColor', cols(s, :), 'MarkerEdgeColor', 'k');
    end
end
xlabel(a, 'horizontal distance from the origin (m)'); ylabel(a, 'Z height (m)'); title(a, 'Side view (height)');
% measured vs true distance
a = nexttile(tl); hold(a, 'on'); grid(a, 'on'); box(a, 'on');
mx = max([R(ok, 4); R(ok, 5); 1]) * 1.15;
plot(a, [0 mx], [0 mx], 'k-', 'LineWidth', 1);
text(a, 0.62 * mx, 0.68 * mx, 'perfect', 'Rotation', 45, 'Color', [0.4 0.4 0.4]);
for s = find(ok)'
    plot(a, R(s, 4), R(s, 5), 'o', 'MarkerSize', 10, 'MarkerFaceColor', cols(s, :), 'MarkerEdgeColor', 'k');
    text(a, R(s, 4), R(s, 5), sprintf('  %d', s), 'Color', cols(s, :), 'FontWeight', 'bold');
end
xlim(a, [0 mx]); ylim(a, [0 mx]); axis(a, 'square');
xlabel(a, 'true distance (m)'); ylabel(a, 'measured distance (m)'); title(a, 'Distance: measured vs true');
% direction error and 3-D error per point
a = nexttile(tl); hold(a, 'on'); grid(a, 'on'); box(a, 'on');
yyaxis(a, 'left');
bar(a, find(ok), R(ok, 6), 0.4, 'FaceColor', [0.2 0.45 0.75]);
ylabel(a, 'direction error (deg)');
yyaxis(a, 'right');
plot(a, find(ok), 100 * vecnorm(R(ok, 1:3) - test(ok, :), 2, 2), 'o-', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto');
ylabel(a, '3-D position error (cm)');
xticks(a, 1:nTest); xlabel(a, 'test point');
title(a, 'Errors per test point (bars: direction, line: 3-D position)');
png = fullfile(outdir, [stamp '_tests.png']);
exportgraphics(fig, png, 'Resolution', 110);
fprintf('picture saved to %s\n', png);
end
