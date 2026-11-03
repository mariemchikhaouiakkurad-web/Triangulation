function tdoa_live(cfg)
% TDOA_LIVE  Real-time 3-D localisation (x, y, z) with the own TDOA method (all 3 boards).
%
%   Every cfg.block s:  tdoa_estimate on every board (142 pairs, no sync between boards)
%                       -> tdoa_locate / tdoa_solve -> position + 95 % region (tdoa_region)
%   Display: tdoa_view (top view out to cfg.tdoa.max_range, zoom, side view, azimuth /
%   elevation / distance over time with 95 % bands).
%   CALIBRATION (tdoa_calibrate, 3-D pose of every board): the phone plays white noise at
%   the points cfg.tdoa.cal.  It is saved to tdoa_calibration.mat and reused on the next
%   run while cfg.tdoa.reuse_cal = true (recalibrate after moving a board).
%   Close the figure to stop.  Called from tri_config.m when cfg.mode = 'tdoa_live'.
if nargin < 1                                        % typed alone: take the settings from tri_config
    fprintf('No settings given - running tri_config (set cfg.mode = ''tdoa_live'' there).\n');
    tri_config;
    return
end
fs = cfg.fs; frame = round(0.1 * fs); blk = round(cfg.block * fs);
G = 1; if isfield(cfg, 'mic_gain'), G = cfg.mic_gain; end    % software amplifier on all mics
S = cfg.tdoa;
B = tdoa_boards(cfg);
here = fileparts(mfilename('fullpath'));
calfile = fullfile(here, 'tdoa_calibration.mat');
cal = S.cal; if size(cal, 2) == 2, cal(:, 3) = S.src_z; end

V = tdoa_view('init', cfg, B);

% ------------------------------------------------------------------ audio + session files
dev = tri_devices();
outdir = fullfile(here, 'live_sessions');
if ~isfolder(outdir), mkdir(outdir); end
stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
rd = cell(1, 3); wr = cell(1, 3);
for b = 1:3
    rd{b} = audioDeviceReader('Driver', dev(b).driver, 'Device', dev(b).name, ...
        'NumChannels', dev(b).nch, 'SampleRate', fs, 'SamplesPerFrame', frame, 'OutputDataType', 'double');
    setup(rd{b});
    wr{b} = dsp.AudioFileWriter(fullfile(outdir, [stamp '_' dev(b).file]), 'SampleRate', fs, 'DataType', 'int32');
end
logf = fopen(fullfile(outdir, [stamp '_tdoa_live.csv']), 'w');
if logf > 2
    fprintf(logf, ['t_s,level_dBFS,x,y,z,std_x,std_y,std_z,rms_us,conf_A,conf_U8,conf_B,r95_min,r95_mid,r95_max,' ...
                   'azimuth_deg,elevation_deg,distance_m\n']);
end
cleanup = onCleanup(@() close_all(rd, wr, logf)); %#ok<NASGU>
fprintf('Saving this session to %s (%s_*)\n', outdir, stamp);
for k = 1:10, for b = 1:3, rd{b}(); end, end         % warm-up
buf = {[], [], []}; t0 = tic;

% ------------------------------------------------------------------ calibration
if S.reuse_cal && isfile(calfile)
    L = load(calfile);
    for b = 1:3, B(b).pose = L.pose(b); end
    fprintf('Calibration loaded from %s (%s). Set cfg.tdoa.reuse_cal = false after moving a board.\n', ...
        calfile, L.when);
    if isfield(L, 'sys')
        sys = L.sys;
    else
        sys = [];
        fprintf('  (old calibration file: no systematic error measured - the 95 %% region shows the noise only)\n');
    end
else
    nCal = size(cal, 1);
    taucal = arrayfun(@(b) zeros(size(B(b).pairs, 1), 0), 1:3, 'UniformOutput', false);
    taucal = repmat(taucal, nCal, 1);                 % {point, board} -> pairs x blocks
    for p = 1:nCal
        tdoa_view('target', V, p);                    % red circle on the point to use now
        tm = tic;
        while toc(tm) < S.move_sec && ishghandle(V.fig)
            [buf, ~] = grab(rd, wr, cfg.use_ch, buf, blk, G);
            buf = {[], [], []};
            tdoa_view('status', V, sprintf('CALIBRATION %d/%d: phone at cal %d (x=%.2f, y=%.2f, z=%.2f m), white noise ON - %.0f s', ...
                p, nCal, p, cal(p, :), S.move_sec - toc(tm)));
        end
        tc = tic;
        while toc(tc) < S.rec_sec && ishghandle(V.fig)
            [buf, X] = grab(rd, wr, cfg.use_ch, buf, blk, G);
            tdoa_view('status', V, sprintf('CALIBRATION %d/%d: keep the phone STILL at cal %d - %.0f s', ...
                p, nCal, p, S.rec_sec - toc(tc)));
            if isempty(X) || toc(tc) < S.settle_sec, continue; end
            lev = 20 * log10(rms(X{1}(:)) + eps);
            if lev < cfg.live_gate, continue; end
            for b = 1:3
                taucal{p, b}(:, end + 1) = tdoa_estimate(X{b} - mean(X{b}, 1), B(b).pairs, B(b).maxlag, fs, S);
            end
        end
        if ~ishghandle(V.fig), return; end
        if isempty(taucal{p, 1})
            error(['cal %d: the phone was not heard (nothing louder than cfg.live_gate = %d dBFS).\n' ...
                   'White noise at FULL volume, speaker facing the boards, and keep it playing.'], p, cfg.live_gate);
        end
        fprintf('cal %d done (%d blocks)\n', p, size(taucal{p, 1}, 2));
    end
    for b = 1:3
        T = cell2mat(cellfun(@(z) median(z, 2), taucal(:, b)', 'UniformOutput', false));
        [B(b).pose, e] = tdoa_calibrate(B(b).Mloc, B(b).centre, B(b).pairs, T, cal, S.src_z, S.flat(b), prior_of(S, b));
        tdoa_print_pose(B(b), e);
    end
    sys = tdoa_systematic(B, cellfun(@(z) median(z, 2), taucal, 'UniformOutput', false), cal, cfg);
    pose = [B.pose]; when = char(datetime('now')); %#ok<NASGU>
    save(calfile, 'pose', 'when', 'sys');
    fprintf('Calibration saved to %s\n', calfile);
end
tdoa_view('target', V, 0);
tdoa_view('mics', V, B);

% ------------------------------------------------------------------ live
fprintf('\n=== LIVE TDOA: play a sound anywhere. Close the figure to stop. ===\n');
o = V.o;
while ishghandle(V.fig) && toc(t0) < cfg.live_seconds
    [buf, X] = grab(rd, wr, cfg.use_ch, buf, blk, G);
    if isempty(X), continue; end
    tnow = toc(t0);
    lev = 20 * log10(rms(X{1}(:)) + eps);
    if lev < cfg.live_gate
        tdoa_view('status', V, sprintf('quiet (%.0f dBFS) - play a sound from the phone', lev));
        continue;
    end
    tau = cell(1, 3); conf = tau;
    for b = 1:3
        [tau{b}, conf{b}] = tdoa_estimate(X{b} - mean(X{b}, 1), B(b).pairs, B(b).maxlag, fs, S);
    end
    [p, info] = tdoa_locate(B, tau, conf, cfg);
    cf = cellfun(@mean, conf);
    good = info.rms_us < S.max_rms_us;
    [Ct, r95] = tdoa_region(p, info.cov, sys);       % 95 % region: noise + systematic error
    d = p - o; az = atan2d(d(1), d(2)); el = asind(d(3) / max(norm(d), eps));
    fprintf(['t=%6.1f  az %+5.0f  el %+4.0f deg  dist %5.2f m  |  x=%+.2f y=%+.2f z=%+.2f m' ...
             '  (95%% region +-%.0f/%.0f/%.0f cm)  fit %4.1f us  [%.0f dBFS]%s\n'], ...
        tnow, az, el, norm(d), p, 100 * r95, info.rms_us, lev, ternary(good, '', '  (uncertain)'));
    if logf > 2
        fprintf(logf, '%.2f,%.1f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.1f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.1f,%.1f,%.3f\n', ...
            tnow, lev, p, info.std_xy, info.rms_us, cf, r95, az, el, norm(d));
    end
    V = tdoa_view('update', V, tnow, p, Ct, good, lev);
end
fprintf('Stopped.\n');
end

% =====================================================================================
function [buf, X] = grab(rd, wr, use_ch, buf, blk, G)
% read 1 frame per board; X = last block of every board once blk new samples arrived
for b = 1:numel(rd)
    [x, lost] = rd{b}();
    tri_write(wr{b}, x, lost);
    buf{b} = [buf{b}; zeros(double(lost), numel(use_ch{b})); x(:, use_ch{b})];
end
X = {};
if all(cellfun(@(z) size(z, 1), buf) >= blk)
    X = cellfun(@(z) G * z(end - blk + 1:end, :), buf, 'UniformOutput', false);   % amplified
    buf = {[], [], []};
end
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
