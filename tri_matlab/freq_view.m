function V = freq_view(action, V, varargin)
% FREQ_VIEW  Figure of the frequency (harmonic comb) detection and its direction.
%
%   V = freq_view('init', cfg, B)
%   freq_view('status', V, text)
%   V = freq_view('update', V, t, det, dir, drone, level)     det from freq_detect, dir from
%                                                             freq_direction ([] when nothing found)
% Panels: spectrum (whitened) with the detected comb | f0 over time | az-el power map
%         top view out to cfg.tdoa.max_range with the direction ray + cone | azimuth | elevation
switch action
    case 'init', V = init(V, varargin{:});
    case 'status', set(V.head, 'String', varargin{1}); drawnow limitrate;
    case 'update', V = update(V, varargin{:});
end
end

% =====================================================================================
function V = init(cfg, B)
F = cfg.freq;
R = cfg.tdoa.max_range; if isempty(R), R = 10; end
P = vertcat(cfg.tri.pos{:}); if size(P, 2) == 2, P(:, 3) = 0; end
V.R = R; V.o = mean(P, 1); V.F = F;
V.red = [0.82 0.15 0.12]; V.grey = [0.6 0.6 0.6];
cols = {[0.12 0.47 0.71], [0.93 0.5 0.05], [0.55 0.35 0.75]};
V.fig = figure('Name', 'Frequency (drone) detection - close to stop', 'Color', 'w', 'Position', [40 40 1500 880]);
set(V.fig, 'DefaultAxesFontSize', 11);
tl = tiledlayout(V.fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
tl.OuterPosition = [0 0 1 0.95];
V.head = annotation(V.fig, 'textbox', [0 0.955 1 0.045], 'String', 'starting ...', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'none');

% 1. spectrum
a = nexttile(tl, 1); hold(a, 'on'); grid(a, 'on'); box(a, 'on');
V.spec = plot(a, nan, nan, '-', 'Color', [0.3 0.3 0.3]);
V.comb = plot(a, nan, nan, 'v', 'MarkerSize', 9, 'MarkerFaceColor', V.red, 'MarkerEdgeColor', 'k');
yline(a, F.peak_db, '--', sprintf('%g dB', F.peak_db), 'Color', V.red, 'LabelHorizontalAlignment', 'left');
xlim(a, [0 F.fmax / 1000]); ylim(a, [-5 40]);
xlabel(a, 'frequency (kHz)'); ylabel(a, 'dB above the noise floor');
title(a, 'Spectrum (all mics), detected harmonics'); V.axs = a;

% 2. f0 over time
a = nexttile(tl, 2); hold(a, 'on'); grid(a, 'on'); box(a, 'on');
V.f0det = plot(a, nan, nan, 'o', 'MarkerSize', 6, 'MarkerFaceColor', V.grey, 'MarkerEdgeColor', 'none');
V.f0dr = plot(a, nan, nan, 'o', 'MarkerSize', 7, 'MarkerFaceColor', V.red, 'MarkerEdgeColor', 'k');
ylim(a, F.f0_range); xlabel(a, 'time (s)'); ylabel(a, 'f0 (Hz)');
title(a, 'Fundamental f0 (red = confirmed drone)'); V.axf = a;

% 3. az-el power map
a = nexttile(tl, 3); hold(a, 'on'); box(a, 'on');
azg = -180:F.az_step:180 - F.az_step; elg = F.el_range(1):F.az_step:F.el_range(2);
V.img = imagesc(a, azg, elg, zeros(numel(elg), numel(azg)));
set(a, 'YDir', 'normal'); clim(a, [0 1]); colormap(a, parula); colorbar(a);
V.mpk = plot(a, nan, nan, 'w+', 'MarkerSize', 16, 'LineWidth', 2.5);
axis(a, 'tight'); xlabel(a, 'azimuth (deg)  0 = toward UMA-8, +90 = right'); ylabel(a, 'elevation (deg)');
title(a, 'Direction map at the detected frequencies');

% 4. top view
a = nexttile(tl, 4); hold(a, 'on'); axis(a, 'equal'); box(a, 'on');
th = linspace(0, 2 * pi, 181);
for r = [1 2 5 10 20 50]
    if r > R, break; end
    plot(a, V.o(1) + r * sin(th), V.o(2) + r * cos(th), ':', 'Color', [0.6 0.6 0.6]);
    text(a, V.o(1) + r * sind(45), V.o(2) + r * cosd(45), sprintf(' %g m', r), 'Color', [0.45 0.45 0.45], 'FontSize', 9);
end
for b = 1:numel(B)
    plot(a, P(b, 1), P(b, 2), 's', 'MarkerSize', 8, 'Color', cols{b}, 'MarkerFaceColor', cols{b});
end
V.cone = patch(a, nan, nan, V.red, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
V.ray = plot(a, nan, nan, '-', 'Color', V.red, 'LineWidth', 2);
xlim(a, V.o(1) + [-R R]); ylim(a, V.o(2) + [-R R]);
xlabel(a, 'X (m)  right'); ylabel(a, 'Y (m)  toward UMA-8'); title(a, 'Direction of the sound (distance unknown)');

% 5-6. azimuth and elevation over time with bands
names = {'Azimuth (deg)', 'Elevation (deg)'}; lims = {[-180 180], F.el_range};
for k = 1:2
    a = nexttile(tl, 4 + k); hold(a, 'on'); grid(a, 'on'); box(a, 'on');
    V.band(k) = patch(a, nan, nan, V.red, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
    V.line(k) = plot(a, nan, nan, '.-', 'Color', V.red, 'MarkerSize', 10);
    ylim(a, lims{k}); xlabel(a, 'time (s)'); title(a, names{k}); V.axt(k) = a;
end
V.hist = zeros(0, 8);                                % t f0 drone az daz el del found
drawnow;
end

function V = update(V, t, det, dir, cls, level)
drone = strcmp(cls, 'DRONE');
set(V.spec, 'XData', det.f / 1000, 'YData', det.snr_db);
set(V.comb, 'XData', det.harm_f / 1000, 'YData', det.harm_snr);
found = ~isempty(dir) && ~isnan(dir.az) && ~strcmp(cls, 'electrical');
if found
    az = dir.az; el = dir.el; daz = dir.daz; del = dir.del;
    set(V.img, 'CData', dir.map); set(V.mpk, 'XData', az, 'YData', el);
    col = V.red; if dir.ambiguous || ~drone, col = V.grey; end
    A = az + linspace(-daz, daz, 40); R = V.R;
    set(V.cone, 'XData', V.o(1) + [0, R * sind(A), 0], 'YData', V.o(2) + [0, R * cosd(A), 0], 'FaceColor', col);
    set(V.ray, 'XData', V.o(1) + [0 R * sind(az)], 'YData', V.o(2) + [0 R * cosd(az)], 'Color', col);
else
    az = NaN; el = NaN; daz = NaN; del = NaN;
end
V.hist(end + 1, :) = [t, det.f0, drone, az, daz, el, del, found];
V.hist = V.hist(max(1, end - 240):end, :);
H = V.hist;
set(V.f0det, 'XData', H(:, 1), 'YData', H(:, 2));
set(V.f0dr, 'XData', H(H(:, 3) > 0, 1), 'YData', H(H(:, 3) > 0, 2));
xlim(V.axf, [max(0, t - 60), t + 1]);
G = H(H(:, 8) > 0, :);
for k = 1:2
    c = 2 + 2 * k;                                    % 4 (az) or 6 (el)
    if ~isempty(G)
        set(V.band(k), 'XData', [G(:, 1)', fliplr(G(:, 1)')], 'YData', [G(:, c)' - G(:, c + 1)', fliplr(G(:, c)' + G(:, c + 1)')]);
        set(V.line(k), 'XData', G(:, 1), 'YData', G(:, c));
    end
    xlim(V.axt(k), [max(0, t - 60), t + 1]);
end
if ~det.detected
    txt = sprintf('no harmonic comb (%s)   [%.0f dBFS]', det.reason, level);
elseif strcmp(cls, 'electrical')
    txt = sprintf('ELECTRICAL tone f0 = %.0f Hz (same phase on every board) - ignored from now on   [%.0f dBFS]', det.f0, level);
elseif ~found
    txt = sprintf('%s  f0 = %.0f Hz (%d harmonics) - direction not determined (%s)   [%.0f dBFS]', ...
        upper(cls), det.f0, numel(det.harm_f), dir.reason, level);
else
    lab = 'candidate'; if drone, lab = 'DRONE'; end
    amb = ''; if dir.ambiguous, amb = '  (ambiguous)'; end
    txt = sprintf('%s  f0 = %.0f Hz (%d harmonics)    azimuth %+.0f%s +-%.0f%s    elevation %+.0f%s +-%.0f%s%s   [%.0f dBFS]', ...
        lab, det.f0, numel(det.harm_f), az, char(176), daz, char(176), el, char(176), del, char(176), amb, level);
end
set(V.head, 'String', txt);
drawnow limitrate;
end
