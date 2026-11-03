function tri_waves_live(cfg)
% TRI_WAVES_LIVE  See the sound waves arrive on every mic of the 3 boards, in real time.
%   cfg.mode = 'waves' in tri_config.   Settings: cfg.wave.
%
% One column per board (UMA16A, UMA-8, UMA16B):
%   row 1  RAW waves, band-passed, one trace per mic, sorted by arrival time.
%          The black dotted line joins the measured arrival of each mic = the WAVEFRONT
%          sweeping across the board.  Its slope is the TDOA.
%   row 2  DELAY-AND-SUM: every trace shifted back by its own delay -> they overlap;
%          thick black = their mean (the beam).  "gain" = how much the sum gains over one mic
%          (perfect alignment: 10*log10(n) dB = 12 dB for 16 mics, 8.5 dB for 7, 3 dB for 2).
%   row 3  DELAY per mic (us) relative to the first mic: bars = measured (GCC from all
%          pairs, least squares), red circles = predicted from the direction found by all
%          3 boards (only with a calibration file).  Bars on circles = calibration + TDOA agree.
%   row 4  PHASE per mic at the strongest frequency f* (what freq mode uses):
%          lines = measured phase vs the first mic, dots = -2*pi*f*tau from the measured
%          delay.  Above ~4 kHz the phase wraps (spatial aliasing): dots and lines still agree
%          mod 360 deg, but the phase alone can no longer tell the delay.
%
% IMPORTANT: the 3 boards have NO common clock.  Inside one column, times are comparable
% to the microsecond; between columns they are not (each column has its own time 0 =
% its own first mic).  The boards are combined through the DIRECTION, not through time.
%
% Keys: SPACE = freeze / resume (the audio keeps flowing), close the figure to stop.
W = cfg.wave; fs = cfg.fs; c = 343;
B = tdoa_boards(cfg); nb = numel(B);

Pw = {};                                          % world mic positions (calibration), optional
try
    pose = tdoa_load_poses(W.cal_file, nb, cfg.tri.pos);
    Pw = cell(1, nb);
    for b = 1:nb, Pw{b} = tdoa_pose(B(b).Mloc, pose(b)); end
catch err
    fprintf('no calibration (%s): predicted delays not shown\n', err.message);
    Pw = {};
end
opt = struct('nfft', W.nfft, 'band', W.band, 'bands', W.band, 'up', 8, 'legacy', false, 'noise', []);
grid = struct('az', -180:4:176, 'el', -10:10:80, 'r', Inf, 'fine_az', 1, 'fine_el', 2);
nblk = round(W.block * fs); nhop = round(W.hop * fs);
half = round(W.win_ms / 2e3 * fs);

lbl = cell(1, nb);
for b = 1:nb
    if b == 3
        lbl{b} = arrayfun(@(k) sprintf('MIC%d', k), cfg.tri.b_mics, 'UniformOutput', false);
    else
        lbl{b} = arrayfun(@(k) sprintf('ch%d', k), B(b).ch, 'UniformOutput', false);
    end
end

src = wave_source(cfg, W, nhop);
cl = onCleanup(@() src.close()); %#ok<NASGU>
[fig, H] = make_fig(B, lbl, W);

buf = cell(1, nb); chs = cell(1, nb);
while isvalid(fig)
    [x, ok] = src.read();
    if ~ok, break; end
    for b = 1:nb
        if isempty(chs{b})
            if size(x{b}, 2) >= max(B(b).ch), chs{b} = B(b).ch; else, chs{b} = 1:numel(B(b).ch); end
            buf{b} = zeros(0, numel(chs{b}));
        end
        buf{b} = [buf{b}; x{b}(:, chs{b})];
        if size(buf{b}, 1) > nblk, buf{b} = buf{b}(end - nblk + 1:end, :); end
    end
    if any(cellfun(@(z) size(z, 1), buf) < nblk), continue; end
    if ~isvalid(fig), break; end
    if getappdata(fig, 'paused'), drawnow; continue; end

    info = cell(1, nb); tm = cell(1, nb);
    for b = 1:nb
        [tau, conf, info{b}] = tdoa_estimate(buf{b}, B(b).pairs, B(b).maxlag, fs, opt);
        tm{b} = mic_delays(tau, conf, B(b).pairs, size(buf{b}, 2));
    end
    u = []; dtxt = '';
    if ~isempty(Pw)
        out = tdoa_srp(info, Pw, B, grid);
        if isfinite(out.az)
            u = [sind(out.az) * cosd(out.el), cosd(out.az) * cosd(out.el), sind(out.el)];
            dtxt = sprintf('   direction from all 3 boards: az %.0f deg, el %.0f deg', out.az, out.el);
        end
    end
    lev = -inf;
    for b = 1:nb
        tp = [];
        if ~isempty(u)
            tp = -(Pw{b} * u') / c; tp = tp - tp(1);   % far-field predicted arrival vs first mic
        end
        lev = max(lev, draw_board(H(b), buf{b}, tm{b}, tp, fs, W, half));
    end
    set(H(1).sg, 'String', sprintf('%.0f dBFS (band %g-%g Hz)%s     [SPACE = freeze]', ...
        lev, W.band(1), W.band(2), dtxt));
    drawnow limitrate
end
end

% ------------------------------------------------------------------------------------------
function t = mic_delays(tau, conf, pairs, n)
% arrival time of every mic vs mic 1 from all pair TDOAs (tau = t_i - t_j), weighted LS
P = size(pairs, 1);
A = zeros(P, n);
A(sub2ind([P n], (1:P)', pairs(:, 1))) = 1;
A(sub2ind([P n], (1:P)', pairs(:, 2))) = -1;
w = sqrt(max(conf(:), 1e-3));
t = [0; (A(:, 2:end) .* w) \ (tau(:) .* w)];
end

function lev = draw_board(h, X, tm, tp, fs, W, half)
[nX, n] = size(X);
f2 = (0:nX - 1)' * fs / nX; f2(f2 > fs / 2) = f2(f2 > fs / 2) - fs;
mask = abs(f2) >= W.band(1) & abs(f2) <= W.band(2);
Y = fft(X) .* mask;
y = real(ifft(Y));                                       % band-passed raw
yal = real(ifft(Y .* exp(1j * 2 * pi * f2 * tm(:)')));    % each mic advanced by its delay
lev = 20 * log10(sqrt(mean(y(:) .^ 2)) + eps);

% window centred on the loudest moment of the aligned sum
m0 = half + ceil(max(abs(tm)) * fs) + 2;
e = movmean(sum(yal, 2) .^ 2, max(3, round(0.002 * fs)));
e([1:min(2 * m0, nX), max(1, nX - 2 * m0 + 1):nX]) = 0;
[~, k] = max(e);
k = min(max(k, 2 * m0 + 1), nX - 2 * m0);
seg = (k - 2 * m0):(k + 2 * m0);                         % wide, then cropped (interpft edges)
U = W.up; L = numel(seg);
ys = interpft(y(seg, :), U * L); ya = interpft(yal(seg, :), U * L);
tt = ((0:U * L - 1)' / (U * fs) - 2 * m0 / fs) * 1e3;    % ms, 0 = window centre
keep = abs(tt) <= half / fs * 1e3;
tt = tt(keep); ys = ys(keep, :); ya = ya(keep, :);
nrm = max(abs(ys), [], 1) + eps;

[~, ord] = sort(tm); rank = zeros(n, 1); rank(ord) = 1:n;
for m = 1:n
    set(h.raw(m), 'XData', tt, 'YData', rank(m) + 0.42 * ys(:, m) / nrm(m));
    set(h.al(m), 'XData', tt, 'YData', ya(:, m) / nrm(m));
end
set(h.wf, 'XData', tm(ord) * 1e3, 'YData', 1:n);
set(h.ax(1), 'YTickLabel', h.lbl(ord));
s = mean(ya ./ nrm, 2);
set(h.sum, 'XData', tt, 'YData', s);
g = 10 * log10(mean(sum(ya, 2) .^ 2) / (sum(mean(ya .^ 2, 1)) + eps) + eps);
g0 = 10 * log10(mean(sum(ys, 2) .^ 2) / (sum(mean(ys .^ 2, 1)) + eps) + eps);
title(h.ax(2), sprintf('delay-and-sum: gain %.1f dB (max %.1f), without delays %.1f dB', ...
    g, 10 * log10(n), g0));

set(h.bar, 'YData', tm(ord) * 1e6);
set(h.ax(3), 'XTickLabel', h.lbl(ord));
if ~isempty(tp)
    set(h.pred, 'XData', 1:n, 'YData', tp(ord) * 1e6);
else
    set(h.pred, 'XData', NaN, 'YData', NaN);
end
lim = max([abs(tm(:)); abs(tp(:)); 1e-5]) * 1e6 * 1.15;
set(h.ax(3), 'YLim', [-lim lim]);

% phase at the strongest frequency in the band
fp = f2 >= W.band(1) & f2 <= W.band(2);
idx = find(fp);
[~, i] = max(mean(abs(Y(idx, :)) .^ 2, 2));
i = idx(i);
ph = angle(Y(i, :) .* conj(Y(i, 1)));
for m = 1:n
    set(h.ph(m), 'XData', [0 cos(ph(m))], 'YData', [0 sin(ph(m))]);
end
pp = -2 * pi * f2(i) * tm(:)';
set(h.pp, 'XData', 1.15 * cos(pp), 'YData', 1.15 * sin(pp));
title(h.ax(4), sprintf('phase vs first mic at f* = %.0f Hz', f2(i)));
end

function [fig, H] = make_fig(B, lbl, W)
nb = numel(B);
fig = figure('Name', 'tri_waves_live', 'NumberTitle', 'off', 'Color', 'w', ...
    'KeyPressFcn', @keyfn);
setappdata(fig, 'paused', false);
tl = tiledlayout(fig, 4, nb, 'TileSpacing', 'compact', 'Padding', 'compact');
sg = title(tl, 'starting...');
H = struct('lbl', {}, 'sg', {}, 'ax', {}, 'raw', {}, 'wf', {}, 'al', {}, 'sum', {}, ...
    'bar', {}, 'pred', {}, 'ph', {}, 'pp', {});
for b = 1:nb
    n = size(B(b).Mloc, 1);
    col = 0.85 * hsv(n);
    h.lbl = lbl{b}; h.sg = sg;
    h.ax = gobjects(1, 4);

    ax = nexttile(tl, b); hold(ax, 'on'); grid(ax, 'on'); h.ax(1) = ax;
    h.raw = gobjects(1, n);
    for m = 1:n, h.raw(m) = plot(ax, NaN, NaN, 'Color', col(m, :), 'LineWidth', 1); end
    h.wf = plot(ax, NaN, NaN, 'k.:', 'MarkerSize', 14, 'LineWidth', 1.2);
    set(ax, 'YLim', [0 n + 1], 'YTick', 1:n, 'YTickLabel', lbl{b}, 'XLim', W.win_ms / 2 * [-1 1]);
    title(ax, sprintf('%s: raw waves (own clock), sorted by arrival', B(b).name));
    xlabel(ax, 'time (ms)');

    ax = nexttile(tl, nb + b); hold(ax, 'on'); grid(ax, 'on'); h.ax(2) = ax;
    h.al = gobjects(1, n);
    for m = 1:n, h.al(m) = plot(ax, NaN, NaN, 'Color', [col(m, :) 0.6], 'LineWidth', 0.8); end
    h.sum = plot(ax, NaN, NaN, 'k', 'LineWidth', 2.5);
    set(ax, 'YLim', [-1.3 1.3], 'XLim', W.win_ms / 2 * [-1 1]);
    xlabel(ax, 'time (ms)');

    ax = nexttile(tl, 2 * nb + b); hold(ax, 'on'); grid(ax, 'on'); h.ax(3) = ax;
    h.bar = bar(ax, 1:n, zeros(1, n), 'FaceColor', [0.55 0.55 0.85]);
    h.pred = plot(ax, NaN, NaN, 'ro', 'MarkerSize', 7, 'LineWidth', 1.5);
    set(ax, 'XTick', 1:n, 'XTickLabel', lbl{b}, 'XLim', [0 n + 1]);
    ylabel(ax, 'delay vs first mic (\mus)');
    title(ax, 'bars = measured, o = predicted from direction');

    ax = nexttile(tl, 3 * nb + b); hold(ax, 'on'); h.ax(4) = ax;
    a = linspace(0, 2 * pi, 100);
    plot(ax, cos(a), sin(a), 'Color', [0.8 0.8 0.8]);
    h.ph = gobjects(1, n);
    for m = 1:n, h.ph(m) = plot(ax, NaN, NaN, 'Color', col(m, :), 'LineWidth', 1.5); end
    h.pp = plot(ax, NaN, NaN, 'k.', 'MarkerSize', 12);
    axis(ax, 'equal'); set(ax, 'XLim', [-1.35 1.35], 'YLim', [-1.35 1.35], 'XTick', [], 'YTick', []);
    box(ax, 'on');
    H(b) = h;
end
end

function keyfn(f, e)
if strcmp(e.Key, 'space')
    p = ~getappdata(f, 'paused');
    setappdata(f, 'paused', p);
    if p, fprintf('frozen (SPACE to resume)\n'); else, fprintf('running\n'); end
end
end

% ------------------------------------------------------------------------------------------
function src = wave_source(cfg, W, nhop)
fs = cfg.fs;
if strcmpi(W.source, 'live')
    dev = tri_devices();
    p = cell(1, numel(dev));
    for b = 1:numel(dev)
        p{b} = audioDeviceReader('Driver', dev(b).driver, 'Device', dev(b).name, ...
            'NumChannels', dev(b).nch, 'SampleRate', fs, 'SamplesPerFrame', nhop, ...
            'OutputDataType', 'double');
        setup(p{b});
    end
    src.read = @() read_live(p);
    src.close = @() close_live(p);
else
    % W.source = a session prefix with a wildcard, e.g. 'live_sessions\20261103_072819_*'
    d = dir(W.source);
    names = lower({d.name});
    key = {'uma16a', 'uma8', 'uma16b'};
    y = cell(1, 3);
    for b = 1:3
        i = find(contains(names, key{b}) & endsWith(names, '.wav'), 1);
        if isempty(i), error('no %s WAV matching %s', key{b}, W.source); end
        [y{b}, fsw] = audioread(fullfile(d(i).folder, d(i).name));
        if fsw ~= fs, error('%s: %d Hz, expected %d', d(i).name, fsw, fs); end
        fprintf('replay %s\n', d(i).name);
    end
    n = min(cellfun(@(z) size(z, 1), y));
    pos = 0;
    src.read = @read_file;
    src.close = @() close_live({});
end
    function [x, ok] = read_file()
        ok = pos + nhop <= n;
        x = {};
        if ~ok, return; end
        x = cellfun(@(z) z(pos + (1:nhop), :), y, 'UniformOutput', false);
        pos = pos + nhop;
        pause(nhop / fs);                         % replay at real speed
    end
end

function close_live(p)
for b = 1:numel(p), release(p{b}); end
end

function [x, ok] = read_live(p)
x = cell(1, numel(p));
for b = 1:numel(p)
    [x{b}, nov] = p{b}();
    if nov > 0, fprintf('  (board %d: %d samples lost)\n', b, nov); end
end
ok = true;
end
