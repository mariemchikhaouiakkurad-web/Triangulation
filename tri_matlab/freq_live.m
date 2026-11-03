function freq_live(cfg, stamp)
% FREQ_LIVE  Drone detection by frequency (harmonic comb) + its direction, all 3 boards.
%
%   Every cfg.freq.hop_sec, on the last cfg.freq.listen_sec of audio:
%     freq_detect    -> fundamental f0 (may drift) + harmonics standing out of the noise
%     freq_direction -> direction (azimuth / elevation) from the harmonics >= cfg.freq.dir_fmin,
%                       and an electrical-tone test per harmonic
%   Classes:  DRONE      comb found in confirm(1) of the last confirm(2) updates, f0 changing by
%                        at most max_f0_jump (relative) between them, not electrical
%             candidate  comb found, not (yet) confirmed
%             electrical most of its harmonics reach all boards with zero delay: their frequencies
%                        are remembered and ignored afterwards (stationary electronic tones)
%   The distance is not estimated (the boards are too close together to measure it far away).
%
%   Calibration: cfg.freq.calibrate = true runs freq_calibrate first (white noise at cfg.tdoa.cal).
%   Otherwise the saved tdoa_calibration.mat is used; it is refused when the board positions in
%   cfg.tri.pos differ from the ones saved with it (a board was moved: recalibrate).
%
%   freq_live(cfg)            live, records live_sessions/<stamp>_*.wav + <stamp>_freq.csv
%   freq_live(cfg, stamp)     the same analysis on a saved session (+ <stamp>_freq.png)
%   Close the figure to stop.  Called from tri_config.m when cfg.mode = 'freq_live'.
if nargin < 1                                        % typed alone: take the settings from tri_config
    fprintf('No settings given - running tri_config (set cfg.mode = ''freq_live'' there).\n');
    tri_config;
    return
end
F = cfg.freq; fs = cfg.fs;
G = 1; if isfield(cfg, 'mic_gain'), G = cfg.mic_gain; end
here = fileparts(mfilename('fullpath'));
outdir = fullfile(here, 'live_sessions');
calfile = fullfile(here, 'tdoa_calibration.mat');
offline = nargin >= 2;
if ~offline && F.calibrate
    freq_calibrate(cfg);
end
if ~isfile(calfile)
    error('No calibration (%s). Set cfg.freq.calibrate = true in tri_config.m and run again.', calfile);
end
L = load(calfile);
if isfield(L, 'pos')
    dpos = max(abs(vertcat(L.pos{:}) - vertcat(cfg.tri.pos{:})), [], 'all');
    if dpos > 0.01 && ~offline
        error(['The board positions in cfg.tri.pos differ from the calibration of %s by %.0f cm: ' ...
               'a board was moved. Set cfg.freq.calibrate = true and run again.'], L.when, 100 * dpos);
    end
else
    fprintf(['NOTE: the calibration of %s has no record of the board positions. If any board was ' ...
             'moved since, set cfg.freq.calibrate = true.\n'], L.when);
end
B = tdoa_boards(cfg);
for b = 1:numel(B), B(b).pose = L.pose(b); end
fprintf('Board orientations from the calibration of %s\n', L.when);
nL = round(F.listen_sec * fs); nH = round(F.hop_sec * fs);
V = freq_view('init', cfg, B);
f0hist = nan(1, F.confirm(2));
blacklist = [];                                      % learned stationary electrical tones (Hz)

if offline
    files = {'uma16a.wav', 'uma8.wav', 'uma16b.wav'};
    x = cell(1, 3);
    for b = 1:3
        y = audioread(fullfile(outdir, [stamp '_' files{b}]), 'native');
        x{b} = G * double(y(:, B(b).ch)) / 2 ^ 31;
    end
    n = min(cellfun(@(z) size(z, 1), x));
else
    stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    dev = tri_devices();
    frame = round(0.1 * fs);
    rd = cell(1, 3); wr = cell(1, 3);
    for b = 1:3
        rd{b} = audioDeviceReader('Driver', dev(b).driver, 'Device', dev(b).name, ...
            'NumChannels', dev(b).nch, 'SampleRate', fs, 'SamplesPerFrame', frame, 'OutputDataType', 'double');
        setup(rd{b});
        wr{b} = dsp.AudioFileWriter(fullfile(outdir, [stamp '_' dev(b).file]), 'SampleRate', fs, 'DataType', 'int32');
    end
end
logf = fopen(fullfile(outdir, [stamp '_freq.csv']), 'w');
if logf > 2
    fprintf(logf, ['t_s,level_dBFS,detected,f0_hz,f0_min_hz,f0_max_hz,n_harm,n_dir_harm,n_elec_harm,score,threshold,' ...
                   'azimuth_deg,az_pm_deg,elevation_deg,el_pm_deg,ambiguity_db,class\n']);
end
if offline
    cleanup = onCleanup(@() fclose_if(logf)); %#ok<NASGU>
    for s = 1:nH:n - nL + 1
        if ~ishghandle(V.fig), break; end
        X = cellfun(@(z) z(s:s + nL - 1, :), x, 'UniformOutput', false);
        [V, f0hist, blacklist] = step(V, X, (s - 1 + nL) / fs, f0hist, blacklist);
    end
    if ishghandle(V.fig)
        exportgraphics(V.fig, fullfile(outdir, [stamp '_freq.png']), 'Resolution', 110);
        fprintf('figure saved to %s\n', fullfile(outdir, [stamp '_freq.png']));
    end
else
    cleanup = onCleanup(@() close_all(rd, wr, logf)); %#ok<NASGU>
    fprintf('Saving this session to %s (%s_*)\n', outdir, stamp);
    for k = 1:10, for b = 1:3, rd{b}(); end, end     % warm-up
    buf = cellfun(@(c) zeros(0, numel(c)), cfg.use_ch, 'UniformOutput', false);
    fresh = 0; t0 = tic;
    fprintf('\n=== FREQUENCY MODE: listening for a harmonic comb. Close the figure to stop. ===\n');
    while ishghandle(V.fig)
        for b = 1:3
            [xb, lost] = rd{b}();
            tri_write(wr{b}, xb, lost);
            buf{b} = [buf{b}; zeros(double(lost), numel(cfg.use_ch{b})); G * xb(:, cfg.use_ch{b})];
            buf{b} = buf{b}(max(1, end - nL + 1):end, :);
        end
        fresh = fresh + frame;
        if fresh < nH || any(cellfun(@(z) size(z, 1), buf) < nL), continue; end
        fresh = 0;
        [V, f0hist, blacklist] = step(V, buf, toc(t0), f0hist, blacklist);
    end
end
fprintf('Stopped.\n');

    function [V, f0hist, blacklist] = step(V, X, t, f0hist, blacklist)
        lev = 20 * log10(rms(X{1}(:)) + eps);
        [det, spec] = freq_detect(X, fs, F, blacklist);
        dir = []; cls = 'none';
        if det.detected
            dir = freq_direction(spec, det, B, F);
            if mean(dir.elec) >= 0.5                  % stationary electronic tone, not a drone
                cls = 'electrical';
                blacklist = unique([blacklist, round(det.harm_f(dir.elec))]);
            else
                cls = 'candidate';
            end
        end
        f0new = NaN; if strcmp(cls, 'candidate'), f0new = det.f0; end
        f0hist = [f0hist(2:end), f0new];
        hits = f0hist(~isnan(f0hist));
        if strcmp(cls, 'candidate') && numel(hits) >= F.confirm(1) ...
                && all(abs(diff(hits)) ./ hits(2:end) <= F.max_f0_jump)
            cls = 'DRONE';
        end
        if det.detected
            fprintf(['t=%6.1f  f0 %6.1f Hz (%4.0f-%4.0f, %2d harm, %d for direction, %d electrical)  ' ...
                     'az %+5.0f +-%-3.0f el %+4.0f +-%-3.0f  %s%s  [%.0f dBFS]\n'], ...
                t, det.f0, min(det.f0_track), max(det.f0_track), numel(det.harm_f), sum(dir.used), sum(dir.elec), ...
                dir.az, dir.daz, dir.el, dir.del, cls, tern(dir.ambiguous && ~isnan(dir.az), ' (ambiguous)', ''), lev);
        end
        if logf > 2
            if isempty(dir), d = [nan(1, 5), 0, 0]; else, d = [dir.az dir.daz dir.el dir.del dir.ambiguity sum(dir.used) sum(dir.elec)]; end
            fprintf(logf, '%.2f,%.1f,%d,%.1f,%.1f,%.1f,%d,%d,%d,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%s\n', ...
                t, lev, det.detected, det.f0, min([det.f0_track NaN]), max([det.f0_track NaN]), numel(det.harm_f), ...
                d(6), d(7), det.score, det.thr, d(1:5), cls);
        end
        V = freq_view('update', V, t, det, dir, cls, lev);
    end
end

% =====================================================================================
function r = tern(c, a, b)
if c, r = a; else, r = b; end
end

function fclose_if(f)
if f > 2, fclose(f); end
end

function close_all(rd, wr, logf)
cellfun(@release, wr);
cellfun(@release, rd);
fclose_if(logf);
fprintf('Session saved in the live_sessions folder\n');
end
