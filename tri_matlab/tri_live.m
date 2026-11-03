function tri_live(cfg, sync)
% TRI_LIVE  Real-time localisation with the 3 arrays, using a PHONE for sync and as test source.
%
%   1. CALIBRATION (cfg.live_cal s): the phone lies at the measured position cfg.spk and
%      plays sync_chirp.wav.  The chirps give clock offset + drift of UMA8 and UMA16B
%      relative to UMA16A (same method as tri_tdoa_combine).
%   2. LIVE: stop the chirp, play any sound from the phone anywhere around the arrays.
%      Every cfg.block s the 3 boards are aligned on the UMA16A timeline and SRP-PHAT
%      (cfg.band, all mics combined) gives azimuth / elevation.  Close the figure to stop.
%
%   Called from tri_config.m when cfg.mode = 'live'.
%   tri_live(cfg, sync) skips the calibration (sync{b} = [drift offset] polynomial, testing).
fs = cfg.fs; c0 = 343;
frame = round(0.1 * fs);
solo = strcmp(cfg.mode, 'live16a');
if solo
    % UMA16A alone: geometry known from the board layout, no sync, no measurements.
    % Directions are in the board frame: az 0 = toward the MIC8-MIC7-MIC10-MIC9 edge,
    % +90 = toward the MIC9-MIC11-MIC13-MIC15 edge, el = angle above the mic face.
    nb = 1;
    M = {tri_board_layout(cfg.odascfg{1}, cfg.use_ch{1})};
    if isfield(cfg, 'front') && ~isempty(cfg.front{1}) && cfg.front{1}(3) < 0
        M{1} = M{1} .* [-1 1 -1];                     % mic face down: turn the board over (about y)
    end
    sync = {[0 0]};
else
    nb = numel(cfg.names);
    [M, N] = tri_geometry(cfg); %#ok<ASGLU>
    if nargin < 2 && any(isnan(cfg.spk))
        error('cfg.spk: measure where the phone lies during calibration (NaN found)');
    end
end

% figure first: opening it takes seconds and would make the streams overflow
fig = figure('Name', 'tri_live - close to stop', 'Position', [80 80 1100 520]);
ax1 = subplot(1, 2, 1); hold(ax1, 'on'); grid(ax1, 'on'); axis(ax1, 'equal');
cols = {[0.12 0.47 0.71], [1 0.5 0.05], [0.58 0.4 0.74]};
for b = 1:nb
    plot(ax1, M{b}(:, 1), M{b}(:, 2), 'o', 'Color', cols{b}, 'MarkerFaceColor', cols{b});
    text(ax1, mean(M{b}(:, 1)), mean(M{b}(:, 2)) - 0.05, cfg.names{b}, 'Color', cols{b}, 'FontWeight', 'bold');
end
ctr = mean(vertcat(M{:}), 1);
xlim(ax1, ctr(1) + [-1.2 1.2]); ylim(ax1, ctr(2) + [-1.2 1.2]);
if solo
    xlim(ax1, [-1.2 1.2]); ylim(ax1, [-1.2 1.2]);
    xlabel(ax1, 'board x (m)  toward MIC9/MIC15 edge'); ylabel(ax1, 'board y (m)  toward MIC8/MIC9 edge');
else
    xlabel(ax1, 'X (m)  right'); ylabel(ax1, 'Y (m)  toward UMA-8');
end
hRay = plot(ax1, ctr(1) * [1 1], ctr(2) * [1 1], 'r-', 'LineWidth', 3);
hTip = plot(ax1, ctr(1), ctr(2), 'ro', 'MarkerSize', 14, 'MarkerFaceColor', 'r');
ax2 = subplot(1, 2, 2); hold(ax2, 'on'); grid(ax2, 'on');
hAz = plot(ax2, nan, nan, 'b.-'); hEl = plot(ax2, nan, nan, 'm.-');
ylim(ax2, [-180 180]); xlabel(ax2, 'time (s)'); ylabel(ax2, 'deg');
legend(ax2, {'azimuth', 'elevation'}, 'Location', 'southwest');
if solo
    title(ax2, 'UMA16A alone: az 0 = MIC8/MIC9 edge, +90 = MIC9/MIC15 edge');
else
    title(ax2, 'az 0 = +Y (toward UMA-8), +90 = +X (right)');
end

dev = tri_devices();
rd = cell(1, nb);
for b = 1:nb
    rd{b} = audioDeviceReader('Driver', dev(b).driver, 'Device', dev(b).name, ...
        'NumChannels', dev(b).nch, 'SampleRate', fs, 'SamplesPerFrame', frame, ...
        'OutputDataType', 'double');
    setup(rd{b});
end
% Save the session (raw audio of the used channels + one CSV row per estimate) for later checks.
% Files go next to this script (MATLAB's current folder may not be writable), one set per run.
outdir = fullfile(fileparts(mfilename('fullpath')), 'live_sessions');
if ~isfolder(outdir), mkdir(outdir); end
stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
wr = cell(1, nb);
for b = 1:nb
    wr{b} = dsp.AudioFileWriter(fullfile(outdir, [stamp '_' dev(b).file]), 'SampleRate', fs, 'DataType', 'int32');
end
logname = fullfile(outdir, [stamp '_log.csv']);
logf = fopen(logname, 'w');
if logf < 0
    warning('Cannot write %s - running without the CSV log', logname);
    logf = 1;                                        % fid 1 = console, so the log lines still work
else
    fprintf(logf, 't_s,az_deg,el_deg,score,level_dBFS\n');
    fprintf('Saving this session to %s\n', outdir);
end
cleanup = onCleanup(@() close_all(rd, wr, logf)); %#ok<NASGU>
for k = 1:10, for b = 1:nb, rd{b}(); end, end       % warm-up: discard start-up backlog

ring = cell(1, nb);
for b = 1:nb
    ring{b} = struct('data', zeros(4 * fs, numel(cfg.use_ch{b})), 'cnt', 0);
end

% ------------------------------------------------------------------ 1. calibration
if nargin < 2 && ~solo
    nCal = ceil(cfg.live_cal * fs / frame);
    syncsig = cell(1, nb);
    for b = 1:nb, syncsig{b} = zeros(nCal * frame + fs, 1); end
    fprintf('\n=== CALIBRATION: phone at cfg.spk playing sync_chirp.wav, %d s ===\n', cfg.live_cal);
    for k = 1:nCal
        for b = 1:nb
            [x, lost] = read_board(rd{b}, wr{b}, cfg.use_ch{b}, dev(b).file);
            n0 = ring{b}.cnt;
            ring{b} = ring_write(ring{b}, x, lost);
            syncsig{b}(n0 + lost + (1:size(x, 1))) = x(:, 1);
        end
        if mod(k, 10) == 0
            fprintf('  %2.0f s left\n', (nCal - k) * frame / fs);
            title(ax1, sprintf('CALIBRATION: phone plays the chirp at cfg.spk - %.0f s left', (nCal - k) * frame / fs));
            drawnow limitrate;
        end
    end
    tpl = sync_chirp_template(fs, cfg.chirp.f0, cfg.chirp.f1, cfg.chirp.dur);
    e = cell(1, nb);
    for b = 1:nb
        arr = detect_chirps(syncsig{b}(1:ring{b}.cnt), tpl, fs, cfg.chirp.period);
        d = norm(cfg.spk - M{b}(1, :));
        e{b} = arr - d / c0 * fs;                     % emission time in this board's samples
        fprintf('%-7s %3d chirps found (phone distance %.2f m)\n', cfg.names{b}, numel(arr), d);
    end
    sync = cell(1, nb); sync{1} = [0 0];
    for b = 2:nb
        [pa, pb] = pair_chirps(e{1}, e{b}, cfg.chirp.period * fs);
        if numel(pa) < 3
            error(['%s: fewer than 3 chirps paired with %s. Phone louder / closer to the arrays, ' ...
                   'and keep it playing during the whole calibration.'], cfg.names{b}, cfg.names{1});
        end
        sync{b} = polyfit(pa, pb - pa, 1);            % board-b index = A index + p(1)*A + p(2)
        rus = std((pb - pa) - polyval(sync{b}, pa)) / fs * 1e6;
        fprintf('%-7s offset %+8.3f ms   drift %+7.2f ppm   fit residual %.1f us (%d chirps)\n', ...
            cfg.names{b}, sync{b}(2) / fs * 1e3, sync{b}(1) * 1e6, rus, numel(pa));
        if rus > 50, fprintf('   WARNING: residual > 50 us: sync unreliable, redo the calibration\n'); end
    end
end

% ------------------------------------------------------------------ 2. live
Mall = vertcat(M{:});
blk = round(cfg.block * fs);
nextEnd = 0; t0 = tic; hist = zeros(0, 3);
if solo
    fprintf('\n=== LIVE (UMA16A alone): play a sound from the phone. Close the figure to stop. ===\n');
else
    fprintf('\n=== LIVE: stop the chirp, play a sound from the phone. Close the figure to stop. ===\n');
end
while ishghandle(fig) && toc(t0) < cfg.live_seconds
    for b = 1:nb
        [x, lost] = read_board(rd{b}, wr{b}, cfg.use_ch{b}, dev(b).file);
        ring{b} = ring_write(ring{b}, x, lost);
    end
    % newest UMA16A block whose aligned samples exist on every board
    nEnd = ring{1}.cnt - 1;
    for b = 2:nb
        p = sync{b};
        nEnd = min(nEnd, floor((ring{b}.cnt - 4 - p(2)) / (1 + p(1))));
    end
    if nEnd < nextEnd, continue; end
    n = (nEnd - blk + 1:nEnd)';
    X = zeros(blk, size(Mall, 1)); col = 0;
    for b = 1:nb
        src = n + polyval(sync{b}, n);                % where A-sample n sits in board b
        i0 = floor(src(1)) - 3; i1 = ceil(src(end)) + 3;
        % if processing fell > 3.5 s behind, samples older than the ring are zero
        seg = ring_read(ring{b}, max(i0, ring{b}.cnt - size(ring{b}.data, 1)), i1);
        seg = [zeros((i1 - i0 + 1) - size(seg, 1), size(seg, 2)); seg];
        nc = size(seg, 2);
        X(:, col + (1:nc)) = interp1((i0:i1)', seg, src, 'spline');
        col = col + nc;
    end
    nextEnd = nEnd + blk;
    lev = 20 * log10(rms(X(:)) + eps);
    tnow = toc(t0);
    if lev < cfg.live_gate
        title(ax1, sprintf('quiet (%.0f dBFS) - play a sound from the phone', lev));
        drawnow limitrate; continue;
    end
    X = X - mean(X, 1);
    [az, el, score] = tri_srp(X, Mall, fs, cfg);
    fprintf('t=%6.1f s   az %+6.1f   el %5.1f deg   score %.2f   level %.0f dBFS\n', tnow, az, el, score, lev);
    if logf > 2, fprintf(logf, '%.2f,%.1f,%.1f,%.3f,%.1f\n', tnow, az, el, score, lev); end
    u = [sind(az) * cosd(el), cosd(az) * cosd(el)];
    set(hRay, 'XData', ctr(1) + [0 u(1)], 'YData', ctr(2) + [0 u(2)]);
    set(hTip, 'XData', ctr(1) + u(1), 'YData', ctr(2) + u(2));
    title(ax1, sprintf('az %+.0f°   el %.0f°   score %.2f', az, el, score));
    hist(end + 1, :) = [tnow az el]; %#ok<AGROW>
    hist = hist(max(1, end - 120):end, :);
    set(hAz, 'XData', hist(:, 1), 'YData', hist(:, 2));
    set(hEl, 'XData', hist(:, 1), 'YData', hist(:, 3));
    xlim(ax2, [max(0, tnow - 30), tnow + 1]);
    drawnow limitrate;
end
fprintf('Stopped.\n');
end

% =====================================================================================
function [x, lost] = read_board(r, w, ch, name)
[x, lost] = r();
lost = double(lost);
x = x(:, ch);
if lost > 0, fprintf('  %s: %d samples dropped (timeline kept with zeros)\n', name, lost); end
tri_write(w, x, lost);                               % saved file keeps the same timeline
end

function close_all(rd, wr, logf)
cellfun(@release, wr);
cellfun(@release, rd);
if logf > 2, fclose(logf); end
fprintf('Session saved in the live_sessions folder\n');
end

function R = ring_write(R, x, lost)
x = [zeros(lost, size(x, 2)); x];                    % keep sample indices continuous
L = size(R.data, 1);
rows = mod(R.cnt + (0:size(x, 1) - 1), L) + 1;
R.data(rows, :) = x;
R.cnt = R.cnt + size(x, 1);
end

function seg = ring_read(R, i0, i1)
seg = R.data(mod(i0:i1, size(R.data, 1)) + 1, :);
end

function arr = detect_chirps(x, tpl, fs, period)
% matched filter -> one arrival (sub-sample) per chirp  (as in tri_tdoa_combine)
N = 2 ^ nextpow2(numel(x) + numel(tpl));
r = real(ifft(fft(x, N) .* conj(fft(tpl, N))));
r = r(1:numel(x));
X = fft(r); h = zeros(size(X)); h(1) = 1; h(2:floor(N / 2)) = 2; h(floor(N / 2) + 1) = 1;
env = abs(ifft(X .* h));
env = env(1:numel(x));
thr = 0.25 * max(env);
gap = round(0.5 * period * fs);
arr = [];
i = 1;
while i <= numel(env)
    if env(i) > thr
        j = min(numel(env), i + gap);
        [~, k] = max(env(i:j)); k = k + i - 1;
        lo = max(2, k - 4); hi = min(numel(r) - 1, k + 4);
        [~, q] = max(r(lo:hi)); q = q + lo - 1;
        a = r(q - 1); b = r(q); c = r(q + 1);
        arr(end + 1, 1) = (q - 1) + 0.5 * (a - c) / (a - 2 * b + c); %#ok<AGROW>
        i = j + 1;
    else
        i = i + 1;
    end
end
end

function [pa, pb] = pair_chirps(ea, eb, per)
[~, j] = min(abs(eb - ea(1)));
c0 = eb(j) - ea(1);
pa = []; pb = [];
for k = 1:numel(ea)
    [d, j] = min(abs(eb - ea(k) - c0));
    if d < 0.1 * per
        pa(end + 1, 1) = ea(k); pb(end + 1, 1) = eb(j); %#ok<AGROW>
    end
end
end
