function V = tdoa_view(action, V, varargin)
% TDOA_VIEW  Figure of the TDOA localisation (used live by tdoa_live and by tdoa_replay).
%
%   V = tdoa_view('init', cfg, B)            create the figure (B = boards, cfg.tdoa.cal shown)
%   tdoa_view('status', V, text)             message in the headline
%   tdoa_view('target', V, p)                highlight calibration point p (0 = none)
%   tdoa_view('mics', V, B)                  draw the calibrated mic positions
%   V = tdoa_view('update', V, t, p, C, good, level)   new estimate p [x y z], 3x3 covariance C
%
% Panels:  top view out to cfg.tdoa.max_range (range rings, direction ray + 95 % cone)
%          top view zoom around the boards       side view (horizontal distance - height)
%          azimuth, elevation and distance over time, each with its 95 % band
% Direction (azimuth, elevation) is what a small array measures well; the distance band
% shows honestly how well the distance is known.
switch action
    case 'init', V = init(V, varargin{:});
    case 'status', set(V.head, 'String', varargin{1}); drawnow limitrate;
    case 'target', target(V, varargin{:});
    case 'mics', mics(V, varargin{:});
    case 'update', V = update(V, varargin{:});
end
end

% =====================================================================================
function V = init(cfg, B)
R = cfg.tdoa.max_range; if isempty(R), R = max(abs([cfg.tri.grid_x cfg.tri.grid_y])); end
P = vertcat(cfg.tri.pos{:}); if size(P, 2) == 2, P(:, 3) = 0; end
cal = cfg.tdoa.cal; if size(cal, 2) == 2, cal(:, 3) = cfg.tdoa.src_z; end
V.R = R; V.o = mean(P, 1); V.cal = cal; V.P = P;
V.cols = {[0.12 0.47 0.71], [0.93 0.5 0.05], [0.55 0.35 0.75]};
V.red = [0.82 0.15 0.12];
V.fig = figure('Name', 'TDOA localisation - close to stop', 'Color', 'w', 'Position', [40 40 1500 880]);
set(V.fig, 'DefaultAxesFontSize', 11, 'DefaultTextFontSize', 11);
tl = tiledlayout(V.fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
V.head = annotation(V.fig, 'textbox', [0 0.955 1 0.045], 'String', 'starting ...', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'none');
tl.OuterPosition = [0 0 1 0.95];
ang = linspace(0, 2 * pi, 181);

% 1. top view out to R
V.ax1 = nexttile(tl, 1); a = V.ax1; hold(a, 'on'); axis(a, 'equal'); box(a, 'on');
for r = ring_list(R)
    plot(a, V.o(1) + r * sin(ang), V.o(2) + r * cos(ang), ':', 'Color', [0.6 0.6 0.6]);
    text(a, V.o(1) + r * sind(45), V.o(2) + r * cosd(45), sprintf(' %g m', r), 'Color', [0.45 0.45 0.45], 'FontSize', 9);
end
board_marks(a, V, B, 1, 2);
plot(a, cal(:, 1), cal(:, 2), 'kx', 'MarkerSize', 9, 'LineWidth', 1.5);
V.cone1 = patch(a, nan, nan, V.red, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
V.ray1 = plot(a, nan, nan, '--', 'Color', V.red, 'LineWidth', 1.5);
V.ell1 = plot(a, nan, nan, '-', 'Color', V.red, 'LineWidth', 1.2);
V.trail1 = plot(a, nan, nan, '.', 'Color', [0.55 0.55 0.55]);
V.dot1 = plot(a, nan, nan, 'o', 'MarkerSize', 11, 'MarkerFaceColor', V.red, 'MarkerEdgeColor', 'k');
xlim(a, V.o(1) + [-R R]); ylim(a, V.o(2) + [-R R]);
xlabel(a, 'X (m)  right'); ylabel(a, 'Y (m)  toward UMA-8'); title(a, sprintf('Top view, out to %g m', R));

% 2. zoom around the boards
V.ax2 = nexttile(tl, 2); a = V.ax2; hold(a, 'on'); axis(a, 'equal'); box(a, 'on'); grid(a, 'on');
Z = 1.5;
board_marks(a, V, B, 1, 2);
plot(a, cal(:, 1), cal(:, 2), 'kx', 'MarkerSize', 10, 'LineWidth', 2);
[lx, ly, lt] = cal_labels(cal);
text(a, lx + 0.04, ly + 0.06, lt, 'FontSize', 9);
V.tgt2 = plot(a, nan, nan, 'o', 'MarkerSize', 30, 'LineWidth', 3, 'Color', V.red);
V.cone2 = patch(a, nan, nan, V.red, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
V.ray2 = plot(a, nan, nan, '--', 'Color', V.red, 'LineWidth', 1.5);
V.ell2 = plot(a, nan, nan, '-', 'Color', V.red, 'LineWidth', 1.2);
V.trail2 = plot(a, nan, nan, '.', 'Color', [0.55 0.55 0.55]);
V.dot2 = plot(a, nan, nan, 'o', 'MarkerSize', 12, 'MarkerFaceColor', V.red, 'MarkerEdgeColor', 'k');
xlim(a, V.o(1) + [-Z Z]); ylim(a, V.o(2) + [-Z Z]);
xlabel(a, 'X (m)'); ylabel(a, 'Y (m)'); title(a, 'Top view, near the boards (zoom)');

% 3. side view: horizontal distance from the boards vs height
V.ax3 = nexttile(tl, 3); a = V.ax3; hold(a, 'on'); box(a, 'on'); grid(a, 'on');
for b = 1:3
    plot(a, hypot(P(b, 1) - V.o(1), P(b, 2) - V.o(2)), P(b, 3), 's', 'MarkerSize', 9, ...
        'Color', V.cols{b}, 'MarkerFaceColor', V.cols{b});
end
plot(a, hypot(cal(:, 1) - V.o(1), cal(:, 2) - V.o(2)), cal(:, 3), 'kx', 'MarkerSize', 9, 'LineWidth', 1.5);
yline(a, 0, 'k-', 'boards', 'LabelHorizontalAlignment', 'left', 'FontSize', 9);
V.tgt3 = plot(a, nan, nan, 'o', 'MarkerSize', 26, 'LineWidth', 3, 'Color', V.red);
V.cone3 = patch(a, nan, nan, V.red, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
V.ray3 = plot(a, nan, nan, '--', 'Color', V.red, 'LineWidth', 1.5);
V.ell3 = plot(a, nan, nan, '-', 'Color', V.red, 'LineWidth', 1.2);
V.dot3 = plot(a, nan, nan, 'o', 'MarkerSize', 11, 'MarkerFaceColor', V.red, 'MarkerEdgeColor', 'k');
xlim(a, [0 R]); ylim(a, [-0.25 * R, 0.75 * R]); daspect(a, [1 1 1]);
xlabel(a, 'horizontal distance from the boards (m)'); ylabel(a, 'Z height (m)');
title(a, 'Side view (height)');

% 4-6. azimuth, elevation, distance over time with 95 % bands
names = {'Azimuth (deg)   0 = toward UMA-8, +90 = right', 'Elevation (deg)   above the boards', ...
         'Distance from the boards (m)'};
lims = {[-180 180], [-30 90], [0 R]};
for k = 1:3
    a = nexttile(tl, 3 + k); hold(a, 'on'); box(a, 'on'); grid(a, 'on');
    V.band(k) = patch(a, nan, nan, V.red, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
    V.line(k) = plot(a, nan, nan, '.-', 'Color', V.red, 'LineWidth', 1.2, 'MarkerSize', 10);
    ylim(a, lims{k}); xlabel(a, 'time (s)'); title(a, names{k}); V.axt(k) = a;
end
V.hist = zeros(0, 10);                               % t az az_lo az_hi el el_lo el_hi r r_lo r_hi
V.trail = zeros(0, 3);
V.ang = ang;
drawnow;
end

function target(V, p)
if p < 1
    set([V.tgt2 V.tgt3], 'XData', nan, 'YData', nan);
else
    c = V.cal(p, :);
    set(V.tgt2, 'XData', c(1), 'YData', c(2));
    set(V.tgt3, 'XData', hypot(c(1) - V.o(1), c(2) - V.o(2)), 'YData', c(3));
end
drawnow limitrate;
end

function mics(V, B)
for b = 1:numel(B)
    W = tdoa_pose(B(b).Mloc, B(b).pose);
    plot(V.ax2, W(:, 1), W(:, 2), '.', 'Color', V.cols{b}, 'MarkerSize', 9);
end
end

function V = update(V, t, p, C, good, level)
d = p - V.o; r = norm(d); e = d / max(r, eps);
az = atan2d(d(1), d(2)); el = asind(max(-1, min(1, d(3) / max(r, eps))));
ea = [cosd(az), -sind(az), 0];                        % horizontal, across the line of sight
ee = cross(e, ea); ee = ee / max(norm(ee), eps);      % vertical, across the line of sight
k = 1.96;                                             % 95 % for one quantity
daz = min(180, k * atan2d(sqrt(max(ea * C * ea', 0)), r));
del = min(90, k * atan2d(sqrt(max(ee * C * ee', 0)), r));
sr = k * sqrt(max(e * C * e', 0));
rlo = max(0, r - sr); rhi = min(V.R, r + sr);
col = V.red; if ~good, col = [0.6 0.6 0.6]; end
R = V.R;

% direction ray + 95 % azimuth cone (top views)
A = az + linspace(-daz, daz, 40);
for a = [V.cone1 V.cone2]
    set(a, 'XData', V.o(1) + [0, R * sind(A), 0], 'YData', V.o(2) + [0, R * cosd(A), 0], 'FaceColor', col);
end
for h = [V.ray1 V.ray2]
    set(h, 'XData', V.o(1) + [0 R * sind(az)], 'YData', V.o(2) + [0 R * cosd(az)], 'Color', col);
end
% 95 % ellipses only while the distance is known; otherwise the cone + distance band show it
showEll = sr < 0.5 * R;
e2 = ell2(C([1 2], [1 2]), V.ang);
if ~showEll, e2(:) = nan; end
set([V.ell1 V.ell2], 'XData', p(1) + e2(1, :), 'YData', p(2) + e2(2, :), 'Color', col);
set([V.dot1 V.dot2], 'XData', p(1), 'YData', p(2), 'MarkerFaceColor', col);
% side view: elevation ray + cone
hdist = hypot(d(1), d(2)); E = el + linspace(-del, del, 30);
set(V.cone3, 'XData', [0, R * cosd(E), 0], 'YData', V.o(3) + [0, R * sind(E), 0], 'FaceColor', col);
set(V.ray3, 'XData', [0 R * cosd(el)], 'YData', V.o(3) + [0 R * sind(el)], 'Color', col);
u = [d(1:2) / max(hdist, eps), 0];                    % horizontal direction of the source
Cs = [u; 0 0 1] * C * [u; 0 0 1]';
e2 = ell2(Cs, V.ang);
if ~showEll, e2(:) = nan; end
set(V.ell3, 'XData', hdist + e2(1, :), 'YData', p(3) + e2(2, :), 'Color', col);
set(V.dot3, 'XData', hdist, 'YData', p(3), 'MarkerFaceColor', col);
if good
    V.trail = [V.trail; p]; V.trail = V.trail(max(1, end - 40):end, :);
    set([V.trail1 V.trail2], 'XData', V.trail(:, 1), 'YData', V.trail(:, 2));
end
% time series with 95 % bands
V.hist(end + 1, :) = [t, az, az - daz, az + daz, el, el - del, el + del, r, rlo, rhi];
V.hist = V.hist(max(1, end - 240):end, :);
H = V.hist; tt = H(:, 1)';
for k = 1:3
    c = 3 * (k - 1) + 2;
    set(V.band(k), 'XData', [tt, fliplr(tt)], 'YData', [H(:, c + 1)', fliplr(H(:, c + 2)')]);
    set(V.line(k), 'XData', tt, 'YData', H(:, c)');
    xlim(V.axt(k), [max(0, t - 60), t + 1]);
end
if sr > 0.5 * R
    rtxt = sprintf('distance %.1f m (95%%: %.1f to >%.0f m - NOT determined)', r, rlo, R);
else
    rtxt = sprintf('distance %.2f m (95%%: %.2f-%.2f m)', r, rlo, rhi);
end
set(V.head, 'String', sprintf(['azimuth %+.0f%s +-%.0f%s    elevation %+.0f%s +-%.0f%s    %s    ' ...
    '|  x %+.2f  y %+.2f  z %+.2f m   [%.0f dBFS]%s'], az, char(176), daz, char(176), el, char(176), del, ...
    char(176), rtxt, p, level, tern(good, '', '  (uncertain)')));
drawnow limitrate;
end

% =====================================================================================
function board_marks(a, V, B, i, j)
for b = 1:3
    plot(a, V.P(b, i), V.P(b, j), 's', 'MarkerSize', 9, 'Color', V.cols{b}, 'MarkerFaceColor', V.cols{b});
    text(a, V.P(b, i), V.P(b, j), ['  ' B(b).name], 'Color', V.cols{b}, 'FontWeight', 'bold', 'FontSize', 9);
end
end

function [x, y, t] = cal_labels(cal)
% one label per spot seen from above: points at the same x, y (only higher) share it,
% e.g. 'cal 1 + cal 4 (40 cm up)'
x = []; y = []; t = {};
for p = 1:size(cal, 1)
    k = find(hypot(x - cal(p, 1), y - cal(p, 2)) < 0.05, 1);
    if isempty(k)
        x(end + 1) = cal(p, 1); y(end + 1) = cal(p, 2); %#ok<AGROW>
        t{end + 1} = sprintf('cal %d', p); %#ok<AGROW>
        if abs(cal(p, 3)) > 0.05, t{end} = sprintf('%s (%.0f cm up)', t{end}, 100 * cal(p, 3)); end
    else
        t{k} = sprintf('%s + cal %d (%.0f cm up)', t{k}, p, 100 * cal(p, 3));
    end
end
end

function r = ring_list(R)
r = [0.5 1 2 5 10 20 50];
r = r(r <= R);
end

function e = ell2(C2, ang)
% 95 % ellipse of a 2 x 2 covariance block (2.45 sigma for two quantities)
[Vv, D] = eig((C2 + C2') / 2);
e = 2.4477 * Vv * sqrt(max(D, 0)) * [cos(ang); sin(ang)];
end

function r = tern(c, a, b)
if c, r = a; else, r = b; end
end
