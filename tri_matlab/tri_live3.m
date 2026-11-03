function tri_live3(cfg)
% TRI_LIVE3  Real-time triangulation with all 3 boards (no clock sync, no chirp needed).
%
%   Each board measures the sound with its own mics only:
%     UMA16A (16 mics), UMA8 (7 mics) -> bearing (SRP-PHAT)
%     UMA16B (MIC1 + MIC16 left)      -> time difference MIC16 - MIC1 (GCC-PHAT, sub-sample)
%   The source position on the table plane is the point that best fits all three
%   (grid search over cfg.tri.grid_x / grid_y).
%
%   1. CALIBRATION: the phone plays white noise at the 2 known points cfg.tri.cal
%      -> rotation + mirror of UMA16A and UMA8, direction of the MIC1->MIC16 axis of UMA16B.
%   2. LIVE: play a sound anywhere. Console + figure + CSV log every cfg.block s.
%   Close the figure to stop.  Called from tri_config.m when cfg.mode = 'live3'.
if nargin < 1                                        % typed alone: take the settings from tri_config
    fprintf('No settings given - running tri_config (set cfg.mode = ''live3'' there).\n');
    tri_config;
    return
end
fs = cfg.fs; c0 = 343;
frame = round(0.1 * fs);
blk = round(cfg.block * fs);
T = cfg.tri;
names = cfg.names;                                    % {'UMA16A', 'UMA8', 'UMA16B'}
Mb = {tri_board_layout(cfg.odascfg{1}, cfg.use_ch{1}), ...
      tri_board_layout(cfg.odascfg{2}, cfg.use_ch{2})};
rB = tri_board_layout(cfg.odascfg{3}, T.b_mics);      % positions of MIC1, MIC16 (read on ch1, ch15)
dB = norm(rB(2, :) - rB(1, :));                       % 42 mm
P = vertcat(T.pos{:}); P = P(:, 1:2);                 % board centres in the world (m), x-y
[gx, gy] = meshgrid(T.grid_x, T.grid_y);
brg = cell(1, 3);                                     % bearing of every grid point from each board
near = false(size(gx));                               % bearings are undefined on top of a board
for b = 1:3
    brg{b} = atan2d(gx - P(b, 1), gy - P(b, 2));
    near = near | hypot(gx - P(b, 1), gy - P(b, 2)) < 0.10;
end

% ------------------------------------------------------------------ figure
fig = figure('Name', 'tri_live3 - close to stop', 'Position', [80 80 1150 560]);
ax1 = subplot(1, 2, 1); hold(ax1, 'on'); grid(ax1, 'on'); axis(ax1, 'equal');
cols = {[0.12 0.47 0.71], [1 0.5 0.05], [0.58 0.4 0.74]};
for b = 1:3
    plot(ax1, P(b, 1), P(b, 2), 's', 'MarkerSize', 12, 'Color', cols{b}, 'MarkerFaceColor', cols{b});
    text(ax1, P(b, 1) + 0.03, P(b, 2) - 0.05, names{b}, 'Color', cols{b}, 'FontWeight', 'bold');
end
plot(ax1, T.cal(:, 1), T.cal(:, 2), 'kx', 'MarkerSize', 12, 'LineWidth', 2);
text(ax1, T.cal(:, 1) + 0.03, T.cal(:, 2), {'cal 1', 'cal 2'});
xlim(ax1, T.grid_x([1 end])); ylim(ax1, T.grid_y([1 end]));
xlabel(ax1, 'X (m)  right'); ylabel(ax1, 'Y (m)  toward UMA-8, away from the laptop');
hRay = gobjects(1, 3);
for b = 1:3, hRay(b) = plot(ax1, nan, nan, '-', 'Color', cols{b}, 'LineWidth', 2); end
hTrail = plot(ax1, nan, nan, '-', 'Color', [0.6 0.6 0.6]);
hSrc = plot(ax1, nan, nan, 'ro', 'MarkerSize', 16, 'MarkerFaceColor', 'r');
ax2 = subplot(1, 2, 2); hold(ax2, 'on'); grid(ax2, 'on');
hAz = [plot(ax2, nan, nan, '.-', 'Color', cols{1}), plot(ax2, nan, nan, '.-', 'Color', cols{2}), ...
       plot(ax2, nan, nan, 'k.-', 'LineWidth', 1.5)];
legend(ax2, {'UMA16A bearing', 'UMA8 bearing', 'source direction'}, 'Location', 'southwest');
ylim(ax2, [-180 180]); xlabel(ax2, 'time (s)'); ylabel(ax2, 'deg (0 = +Y, +90 = +X)');
drawnow;

% ------------------------------------------------------------------ audio + session files
dev = tri_devices();
rd = cell(1, 3); wr = cell(1, 3);
outdir = fullfile(fileparts(mfilename('fullpath')), 'live_sessions');
if ~isfolder(outdir), mkdir(outdir); end
stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
for b = 1:3
    rd{b} = audioDeviceReader('Driver', dev(b).driver, 'Device', dev(b).name, ...
        'NumChannels', dev(b).nch, 'SampleRate', fs, 'SamplesPerFrame', frame, 'OutputDataType', 'double');
    setup(rd{b});
    wr{b} = dsp.AudioFileWriter(fullfile(outdir, [stamp '_' dev(b).file]), 'SampleRate', fs, 'DataType', 'int32');
end
logname = fullfile(outdir, [stamp '_tri_log.csv']);
logf = fopen(logname, 'w');
if logf < 0, warning('Cannot write %s - running without the CSV log', logname); end
if logf > 2
    fprintf(logf, 't_s,phase,level_dBFS,A_az,A_score,U8_az,U8_score,B_tdoa_us,B_angle,src_x,src_y,src_dir,src_dist\n');
end
cleanup = onCleanup(@() close_all(rd, wr, logf)); %#ok<NASGU>
fprintf('Saving this session to %s\n', outdir);
for k = 1:10, for b = 1:3, rd{b}(); end, end         % warm-up

buf = {[], [], []}; t0 = tic;

% ------------------------------------------------------------------ 1. calibration
azcal = cell(2, 2);                                   % {point, board} -> measured board azimuths
taucal = cell(2, 1);                                  % {point} -> UMA16B time differences (s)
for p = 1:2
    t_move = tic;
    while toc(t_move) < T.move_sec && ishghandle(fig)
        title(ax1, sprintf('CALIBRATION %d/2: phone at cal %d (x=%.2f, y=%.2f m), white noise ON - starts in %.0f s', ...
            p, p, T.cal(p, :), T.move_sec - toc(t_move)));
        drawnow limitrate;
        [buf, ~] = grab(rd, wr, cfg.use_ch, buf, blk);  % keep the streams flowing
        buf = {[], [], []};
    end
    fprintf('Calibration point %d: phone at (%.2f, %.2f) m ...\n', p, T.cal(p, :));
    t_cal = tic;
    while toc(t_cal) < T.cal_sec && ishghandle(fig)
        [buf, X] = grab(rd, wr, cfg.use_ch, buf, blk);
        if isempty(X), continue; end
        title(ax1, sprintf('CALIBRATION %d/2: keep the phone STILL at cal %d - %.0f s left', p, p, T.cal_sec - toc(t_cal)));
        drawnow limitrate;
        lev = 20 * log10(rms(X{1}(:)) + eps);
        if lev < cfg.live_gate || toc(t_cal) < T.settle_sec, continue; end   % skip while the phone settles
        for b = 1:2
            azcal{p, b}(end + 1) = tri_srp(X{b} - mean(X{b}, 1), Mb{b}, fs, cfg);
        end
        taucal{p}(end + 1) = gcc_tdoa(X{3}(:, 1), X{3}(:, 2), fs, dB, cfg.band);
        if logf > 2
            fprintf(logf, '%.2f,cal%d,%.1f,%.1f,,%.1f,,%.1f,,,,,\n', toc(t0), p, lev, ...
                azcal{p, 1}(end), azcal{p, 2}(end), taucal{p}(end) * 1e6);
        end
    end
end
if ~ishghandle(fig), return; end

yaw = zeros(1, 2); mir = ones(1, 2);
for b = 1:2
    if isempty(azcal{1, b}) || isempty(azcal{2, b})
        error(['%s heard nothing during calibration (level < cfg.live_gate = %d dBFS). ' ...
               'Phone louder / closer, white noise at full volume.'], names{b}, cfg.live_gate);
    end
    meas = [circmean(azcal{1, b}), circmean(azcal{2, b})];
    truth = [atan2d(T.cal(1, 1) - P(b, 1), T.cal(1, 2) - P(b, 2)), ...
             atan2d(T.cal(2, 1) - P(b, 1), T.cal(2, 2) - P(b, 2))];
    best = inf;
    for m = [1 -1]                                    % -1 = board seen mirrored (mic side down)
        y = circmean(truth - m * meas);
        err = max(abs(wrap180(truth - (m * meas + y))));
        if err < best, best = err; yaw(b) = y; mir(b) = m; end
    end
    spread = max(circstd(azcal{1, b}), circstd(azcal{2, b}));
    fprintf('%-7s rotation %+6.1f deg  %s  calibration error %.1f deg  (estimate spread %.1f deg)\n', ...
        names{b}, yaw(b), ternary(mir(b) < 0, 'mirrored', 'normal  '), best, spread);
    if best > 15 || spread > 15
        fprintf('   WARNING: %s calibration poor - louder noise, phone still, check cfg.tri.pos / cfg.tri.cal\n', names{b});
    end
end

% UMA16B: find the world direction phiB of the MIC1 -> MIC16 axis.
% Plane wave from bearing a:  tdoa(MIC16 - MIC1) = -(dB / c) * cos(a - phiB)
wB = T.b_weight;
if isempty(taucal{1}) || isempty(taucal{2})
    fprintf('UMA16B  heard nothing during calibration - not used\n');
    wB = 0; phiB = 0;
else
    tmeas = [median(taucal{1}), median(taucal{2})];
    truthB = [atan2d(T.cal(1, 1) - P(3, 1), T.cal(1, 2) - P(3, 2)), ...
              atan2d(T.cal(2, 1) - P(3, 1), T.cal(2, 2) - P(3, 2))];
    phis = 0:0.5:359.5;
    pred = -(dB / c0) * cosd(truthB(:) - phis);
    err = sqrt(mean((pred - tmeas(:)) .^ 2, 1));
    [e, i] = min(err); phiB = phis(i);
    spreadB = max(std(taucal{1}), std(taucal{2}));
    fprintf('UMA16B  MIC1->MIC16 axis at %+6.1f deg  calibration error %.0f us  (spread %.0f us, max possible %.0f us)\n', ...
        wrap180(phiB), e * 1e6, spreadB * 1e6, dB / c0 * 1e6);
    if e > 25e-6 || spreadB > 25e-6
        fprintf('   WARNING: UMA16B calibration poor - its weight is reduced\n');
        wB = wB / 4;
    end
end

% ------------------------------------------------------------------ 2. live
fprintf('\n=== LIVE: play a sound anywhere. Close the figure to stop. ===\n');
hist = zeros(0, 4); trail = zeros(0, 2);
while ishghandle(fig) && toc(t0) < cfg.live_seconds
    [buf, X] = grab(rd, wr, cfg.use_ch, buf, blk);
    if isempty(X), continue; end
    tnow = toc(t0);
    lev = 20 * log10(rms(X{1}(:)) + eps);
    if lev < cfg.live_gate
        title(ax1, sprintf('quiet (%.0f dBFS) - play a sound from the phone', lev));
        drawnow limitrate; continue;
    end
    az = zeros(1, 2); sc = zeros(1, 2);
    cost = zeros(size(gx));
    for b = 1:2
        [a, ~, sc(b)] = tri_srp(X{b} - mean(X{b}, 1), Mb{b}, fs, cfg);
        az(b) = wrap180(mir(b) * a + yaw(b));         % world bearing
        cost = cost + sc(b) * (wrap180(brg{b} - az(b)) / T.sigma_deg) .^ 2;
        set(hRay(b), 'XData', P(b, 1) + [0 3 * sind(az(b))], 'YData', P(b, 2) + [0 3 * cosd(az(b))]);
    end
    tauB = gcc_tdoa(X{3}(:, 1), X{3}(:, 2), fs, dB, cfg.band);
    angB = wrap180(phiB + acosd(max(-1, min(1, -tauB * c0 / dB))));   % one of its 2 mirror bearings
    cost = cost + wB * ((tauB + dB / c0 * cosd(brg{3} - phiB)) * 1e6 / T.sigma_us) .^ 2;
    set(hRay(3), 'XData', P(3, 1) + [0 3 * sind(angB)], 'YData', P(3, 2) + [0 3 * cosd(angB)]);

    cost(near) = inf;
    [~, i] = min(cost(:));
    src = [gx(i) gy(i)];
    ctr = mean(P, 1);
    dirs = atan2d(src(1) - ctr(1), src(2) - ctr(2)); dist = norm(src - ctr);
    where = sprintf('source at x=%+.2f y=%+.2f m  (direction %+.0f deg, %.2f m)', src, dirs, dist);
    set(hSrc, 'XData', src(1), 'YData', src(2));
    trail = [trail; src]; trail = trail(max(1, end - 20):end, :); %#ok<AGROW>
    set(hTrail, 'XData', trail(:, 1), 'YData', trail(:, 2));
    fprintf('t=%6.1f  UMA16A %+5.0f (%.2f)  UMA8 %+5.0f (%.2f)  UMA16B %+4.0f us  ->  %s  [%.0f dBFS]\n', ...
        tnow, az(1), sc(1), az(2), sc(2), tauB * 1e6, where, lev);
    if logf > 2
        fprintf(logf, '%.2f,live,%.1f,%.1f,%.3f,%.1f,%.3f,%.1f,%.1f,%.3f,%.3f,%.1f,%.3f\n', ...
            tnow, lev, az(1), sc(1), az(2), sc(2), tauB * 1e6, angB, src, dirs, dist);
    end
    title(ax1, where);
    hist(end + 1, :) = [tnow az dirs]; hist = hist(max(1, end - 120):end, :); %#ok<AGROW>
    for j = 1:3, set(hAz(j), 'XData', hist(:, 1), 'YData', hist(:, j + 1)); end
    xlim(ax2, [max(0, tnow - 30), tnow + 1]);
    drawnow limitrate;
end
fprintf('Stopped.\n');
end

% =====================================================================================
function [buf, X] = grab(rd, wr, use_ch, buf, blk)
% read 1 frame per board; X = last block of every board once blk new samples arrived
for b = 1:numel(rd)
    [x, lost] = rd{b}();
    tri_write(wr{b}, x, lost);
    buf{b} = [buf{b}; zeros(double(lost), numel(use_ch{b})); x(:, use_ch{b})];
end
X = {};
if all(cellfun(@(z) size(z, 1), buf) >= blk)
    X = cellfun(@(z) z(end - blk + 1:end, :), buf, 'UniformOutput', false);
    buf = {[], [], []};
end
end

function tau = gcc_tdoa(x1, x2, fs, d, band)
% time difference x2 - x1 (s) by GCC-PHAT, 16x upsampled for sub-sample resolution
up = 16;
n = 2 ^ nextpow2(2 * numel(x1));
X1 = fft(x1 - mean(x1), n); X2 = fft(x2 - mean(x2), n);
f = (0:n - 1)' * fs / n; f(f > fs / 2) = f(f > fs / 2) - fs;
G = X2 .* conj(X1); G = G ./ (abs(G) + 1e-12);
G(abs(f) < band(1) | abs(f) > band(2)) = 0;
Gu = [G(1:n / 2); zeros((up - 1) * n, 1); G(n / 2 + 1:end)];   % zero-pad the spectrum = upsample
r = real(ifft(Gu));
L = ceil(d / 343 * fs * up) + up;
lags = -L:L;
[~, i] = max(r(mod(lags, n * up) + 1));
tau = lags(i) / (fs * up);
end

function m = circmean(a, w)
if nargin < 2, w = ones(size(a)); end
m = atan2d(sum(w .* sind(a)), sum(w .* cosd(a)));
end

function s = circstd(a)
R = hypot(mean(sind(a)), mean(cosd(a)));
s = sqrt(-2 * log(max(R, 1e-6))) * 180 / pi;
end

function a = wrap180(a)
a = mod(a + 180, 360) - 180;
end

function r = ternary(c, a, b)
if c, r = a; else, r = b; end
end

function close_all(rd, wr, logf)
cellfun(@release, wr);
cellfun(@release, rd);
if logf > 2, fclose(logf); end
fprintf('Session saved in the live_sessions folder\n');
end
