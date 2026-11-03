function lg = tdoa_detect_live(cfg)
% TDOA_DETECT_LIVE  Broadband sound DETECTION + direction with the freq-mode strategy:
%   learn the noise -> only what stands above it counts -> all pairs vote -> CFAR threshold
%   on the vote contrast -> confirm m of n updates with a steady direction.
%
%   cfg.mode = 'tdoa_detect' in tri_config.   Settings: cfg.det (see tri_config).
%   cfg.det.source = 'live'          the 3 boards (tri_devices)
%                  = '<folder>'      offline: uma16a.wav, uma8.wav, uma16b.wav (cfg.files) in it
%
%   Uses the board poses from cfg.det.cal_file (tdoa_calibration.mat), like freq_live.
%   The first cfg.det.learn_sec seconds learn the noise: NO test sound during them.
%   Returns / saves (cfg.det.log_file) one row per update: t, z, thr, det, conf, az, el, r, snr.
D = cfg.det; fs = cfg.fs;
B = tdoa_boards(cfg);
nb = numel(B);
pose = tdoa_load_poses(D.cal_file, nb, cfg.tri.pos);
W = cell(1, nb);
for b = 1:nb, W{b} = tdoa_pose(B(b).Mloc, pose(b)); end

opt = struct('nfft', cfg.tdoa.nfft, 'band', cfg.tdoa.band, 'bands', cfg.tdoa.bands, ...
    'up', cfg.tdoa.up, 'legacy', false, 'noise', []);
grid = struct('az', -180:D.az_step:180 - D.az_step, 'el', cfg.tdoa.el_grid, 'r', D.ranges, ...
    'fine_az', 0.5, 'fine_el', 1);
nblk = round(D.block * fs); nhop = round(D.hop * fs);
nlearn = ceil(D.learn_sec / D.hop);

src = open_source(cfg, D, nhop);
cleanup = onCleanup(@() src.close()); %#ok<NASGU>

buf = cell(1, nb); chs = cell(1, nb);
noise = cell(1, nb);
hist = [];
recent = zeros(0, 4);                                 % [det az el t] of the last updates
lg = struct('t', [], 'z', [], 'thr', [], 'det', [], 'conf', [], 'az', [], 'el', [], 'r', [], 'snr', []);
fig = [];
if D.plot, [fig, hImg, hPt, hTt] = make_plot(grid); end

fprintf('learning the noise for %.0f s: keep the area quiet\n', D.learn_sec);
t = 0; kb = 0;
while t < cfg.live_seconds
    [x, ok] = src.read();
    if ~ok, break; end
    t = t + D.hop;
    for b = 1:nb
        if isempty(chs{b})                            % live: device channels; offline WAV may hold only used ones
            if size(x{b}, 2) >= max(B(b).ch), chs{b} = B(b).ch; else, chs{b} = 1:numel(B(b).ch); end
            buf{b} = zeros(0, numel(chs{b}));
        end
        buf{b} = [buf{b}; x{b}(:, chs{b})];
        if size(buf{b}, 1) > nblk, buf{b} = buf{b}(end - nblk + 1:end, :); end
    end
    if any(cellfun(@(z) size(z, 1), buf) < nblk), continue; end
    kb = kb + 1;

    info = cell(1, nb); snr = nan(1, nb);
    for b = 1:nb
        o = opt; o.noise = noise{b};
        [~, ~, info{b}] = tdoa_estimate(buf{b}, B(b).pairs, B(b).maxlag, fs, o);
        snr(b) = info{b}.snr_db;
    end
    learning = kb <= nlearn;
    haveNoise = ~isempty(noise{1}) && noise{1}.n >= 3;

    out = struct('az', NaN, 'el', NaN, 'r', NaN, 'z', NaN, 'map', []);
    if haveNoise, out = tdoa_srp(info, W, B, grid); end

    % CFAR threshold from the contrast of noise-only blocks (as cfg.freq.cfar_k)
    if numel(hist) >= 8
        m = median(hist); s = 1.4826 * median(abs(hist - m));
        thr = max(D.min_z, m + D.cfar_k * s);
    else
        thr = D.min_z;
    end
    det = ~learning && haveNoise && out.z > thr;

    if ~det                                           % noise-only block: update floor + CFAR history
        for b = 1:nb, noise{b} = tdoa_noise(noise{b}, info{b}, D.noise_alpha); end
        if haveNoise
            hist(end + 1) = out.z; %#ok<AGROW>
            if numel(hist) > D.hist_len, hist = hist(end - D.hist_len + 1:end); end
        end
    end

    % confirm: det in m of the last n updates, all within max_jump deg of the current direction
    recent(end + 1, :) = [det, out.az, out.el, t]; %#ok<AGROW>
    if size(recent, 1) > D.confirm(2), recent = recent(end - D.confirm(2) + 1:end, :); end
    conf = false;
    if det
        dd = recent(recent(:, 1) == 1, :);
        u0 = dir1(out.az, out.el);
        ang = arrayfun(@(k) acosd(min(1, max(-1, dir1(dd(k, 2), dd(k, 3)) * u0'))), 1:size(dd, 1));
        conf = sum(ang <= D.max_jump) >= D.confirm(1);
    end

    if learning, state = 'LEARN';
    elseif conf, state = 'SOUND (confirmed)';
    elseif det, state = 'sound?';
    else, state = '-';
    end
    fprintf('%7.2f s  z=%5.1f thr=%5.1f  SNR=%+5.1f dB  az=%7.1f el=%5.1f r=%4.1f  %s\n', ...
        t, out.z, thr, mean(snr, 'omitnan'), out.az, out.el, out.r, state);

    lg.t(end + 1) = t; lg.z(end + 1) = out.z; lg.thr(end + 1) = thr; lg.det(end + 1) = det;
    lg.conf(end + 1) = conf; lg.az(end + 1) = out.az; lg.el(end + 1) = out.el;
    lg.r(end + 1) = out.r; lg.snr(end + 1) = mean(snr, 'omitnan');

    if D.plot
        if ~isvalid(fig), break; end
        if ~isempty(out.map)
            set(hImg, 'CData', out.map');
            set(hPt, 'XData', out.az, 'YData', out.el);
        end
        set(hTt, 'String', sprintf('%s   z = %.1f (thr %.1f)   az %.1f  el %.1f', state, out.z, thr, out.az, out.el));
        drawnow limitrate
    end
end
if ~isempty(D.log_file)
    save(D.log_file, 'lg', 'cfg');
    fprintf('saved %s\n', D.log_file);
end
end

function u = dir1(az, el)
u = [sind(az) * cosd(el), cosd(az) * cosd(el), sind(el)];
end

function [fig, hImg, hPt, hTt] = make_plot(grid)
fig = figure('Name', 'tdoa_detect', 'NumberTitle', 'off');
hImg = imagesc(grid.az, grid.el, zeros(numel(grid.el), numel(grid.az)));
axis xy; colorbar; hold on;
hPt = plot(NaN, NaN, 'w+', 'MarkerSize', 14, 'LineWidth', 2);
xlabel('azimuth (deg, 0 = +Y ahead, +90 = right)'); ylabel('elevation (deg)');
hTt = title('learning noise');
end

function src = open_source(cfg, D, nhop)
fs = cfg.fs;
if strcmpi(D.source, 'live')
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
    y = cell(1, numel(cfg.files));
    for b = 1:numel(cfg.files)
        [y{b}, fsw] = audioread(fullfile(D.source, cfg.files{b}));
        if fsw ~= fs, error('%s: %d Hz, expected %d', cfg.files{b}, fsw, fs); end
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
