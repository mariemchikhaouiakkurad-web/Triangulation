function res = tri_tdoa_combine(cfg)
% TRI_TDOA_COMBINE  Combine the raw data of 3 arrays (separate USB clocks) into ONE array.
%
%   1. reads the 3 recordings (UMA16A, UMA8, UMA16B)
%   2. finds the sync chirps (6-9 kHz, every 1 s) on each board
%   3. fits clock offset + drift of UMA8 and UMA16B relative to UMA16A, resamples them
%      onto the UMA16A timeline  -> all channels now share one time axis
%   4. builds the world coordinates of every mic (ODAS cfg geometry + board pose)
%   5. SRP-PHAT localisation, 200-4000 Hz, one direction per block:
%         UMA16A alone   vs   all boards combined (cross-board TDOA)
%   6. optional: error against known directions (cfg.truth) + plots
%
%   Run tri_config.m (edit it first).  Works in MATLAB (no toolbox) and GNU Octave.

c0 = 343;
fs = cfg.fs;
nb = numel(cfg.files);
fprintf('\n=== 1. loading recordings ===\n');
x = cell(1, nb); M = cell(1, nb); N = cell(1, nb);
for b = 1:nb
    [y, fsr] = audioread(cfg.files{b});
    if fsr ~= fs, error('%s: fs %d, expected %d', cfg.files{b}, fsr, fs); end
    y = y(:, cfg.use_ch{b});
    x{b} = y - mean(y, 1);
    loc = odas_geometry(cfg.odascfg{b}, cfg.use_ch{b});
    fprintf('%-7s %s: %d ch, %.1f s\n', cfg.names{b}, cfg.files{b}, size(x{b}, 2), size(x{b}, 1) / fs);
    if isfield(cfg, 'geom') && numel(cfg.geom) >= b && ~isempty(cfg.geom{b})
        G = cfg.geom{b};                             % measured mic holes (m) in the common frame
        if any(isnan(G.xyz(:))), error('%s: fill cfg.geom{%d}.xyz with your measurements (NaN found)', cfg.names{b}, b); end
        if size(G.xyz, 1) >= 3
            Lm = odas_geometry(cfg.odascfg{b}, G.ch);
            [Rb, tb, rms_mm] = kabsch(Lm, G.xyz);
            M{b} = (Rb * loc')' + tb;
            N{b} = Rb(:, 3)';                        % board normal; sign set by cfg.front (listening side)
            fr = [0 0 1];
            if isfield(cfg, 'front') && numel(cfg.front) >= b && ~isempty(cfg.front{b}), fr = cfg.front{b}; end
            if dot(N{b}, fr) < 0, N{b} = -N{b}; end
            fprintf('        pose fitted from %d measured mics, fit error %.1f mm%s\n', size(G.xyz, 1), rms_mm, ...
                    ternary(rms_mm > 5, '  <-- > 5 mm: re-check the measurements', ''));
        else                                         % 1-2 mics: use the measured positions directly
            M{b} = zeros(numel(cfg.use_ch{b}), 3);
            for i = 1:numel(cfg.use_ch{b})
                r = find(G.ch == cfg.use_ch{b}(i), 1);
                if isempty(r), error('%s: channel %d has no measured position', cfg.names{b}, cfg.use_ch{b}(i)); end
                M{b}(i, :) = G.xyz(r, :);
            end
            N{b} = [0 0 1];
            fprintf('        %d measured mic positions used directly (spacing %.1f mm)\n', size(M{b}, 1), ...
                    1000 * norm(M{b}(end, :) - M{b}(1, :)));
        end
    else
        Rb = board_rotation(cfg.facing(b), cfg.roll(b), cfg.mount{b});
        M{b} = (Rb * loc')' + cfg.pos{b};
        N{b} = Rb(:, 3)';                            % board normal (front side) in world
    end
end

if isfield(cfg, 'plot') && cfg.plot
    figure('Name', 'geometry check (top view)'); hold on; grid on; axis equal;
    cols = {[0.12 0.47 0.71], [1 0.5 0.05], [0.58 0.4 0.74]};
    for b = 1:nb
        plot(M{b}(:, 1) * 100, M{b}(:, 2) * 100, 'o', 'Color', cols{b}, 'MarkerFaceColor', cols{b});
        text(mean(M{b}(:, 1)) * 100, mean(M{b}(:, 2)) * 100 - 2, cfg.names{b}, 'Color', cols{b}, 'FontWeight', 'bold');
    end
    plot(cfg.spk(1) * 100, cfg.spk(2) * 100, 'kp', 'MarkerSize', 12); text(cfg.spk(1) * 100, cfg.spk(2) * 100 - 2, 'sync speaker');
    xlabel('X (cm)  right'); ylabel('Y (cm)  toward UMA-8'); title('mic positions used - must match the real layout');
end

% ---------------------------------------------------------------- sync
if any(isnan(cfg.spk)), error('cfg.spk: measure the sync speaker position (NaN found)'); end
fprintf('\n=== 2-3. sync chirps -> clock offset + drift ===\n');
tpl = sync_chirp_template(fs, cfg.chirp.f0, cfg.chirp.f1, cfg.chirp.dur);
e = cell(1, nb);
for b = 1:nb
    arr = detect_chirps(x{b}(:, 1), tpl, fs, cfg.chirp.period);
    d = norm(cfg.spk - M{b}(1, :));                 % sync speaker -> first mic of this board
    e{b} = arr - d / c0 * fs;                       % emission time, in this board's samples
    fprintf('%-7s %3d chirps found (speaker distance %.2f m)\n', cfg.names{b}, numel(arr), d);
end
xs = cell(1, nb); xs{1} = x{1};
LA = size(x{1}, 1);
res.sync = struct('name', {}, 'offset_ms', {}, 'drift_ppm', {}, 'resid_us', {});
for b = 2:nb
    [pa, pb] = pair_chirps(e{1}, e{b}, cfg.chirp.period * fs);
    if numel(pa) < 3, error('%s: fewer than 3 chirps paired with %s', cfg.names{b}, cfg.names{1}); end
    p = polyfit(pa, pb - pa, 1);                    % board-b index = A index + p(1)*A + p(2)
    r = (pb - pa) - polyval(p, pa);
    rus = std(r) / fs * 1e6;
    fprintf('%-7s offset %+8.3f ms   drift %+7.2f ppm   fit residual %.1f us (%d chirps)\n', ...
            cfg.names{b}, p(2) / fs * 1e3, p(1) * 1e6, rus, numel(pa));
    if rus > 50
        fprintf('   WARNING: residual > 50 us (~2 samples): check speaker position / chirp level\n');
    end
    res.sync(end + 1) = struct('name', cfg.names{b}, 'offset_ms', p(2) / fs * 1e3, ...
                               'drift_ppm', p(1) * 1e6, 'resid_us', rus);
    n = (0:LA - 1)';
    src = n + polyval(p, n);                        % where A-sample n sits in board b
    xb = x{b};
    xs{b} = zeros(LA, size(xb, 2));
    for m = 1:size(xb, 2)
        xs{b}(:, m) = interp1((0:size(xb, 1) - 1)', xb(:, m), src, 'spline', 0);
    end
end

% ---------------------------------------------------------------- optional yaw calibration
if isfield(cfg, 'calib') && ~isempty(cfg.calib)
    fprintf('\n=== calibration window %.1f-%.1f s, phone at (%.2f %.2f %.2f) m ===\n', cfg.calib(1:2), cfg.calib(3:5));
    i0 = round(cfg.calib(1) * fs) + 1; i1 = min(round(cfg.calib(2) * fs), LA);
    for b = 1:nb
        if size(M{b}, 1) < 3, fprintf('%-7s %d mics: no own direction, skipped\n', cfg.names{b}, size(M{b}, 1)); continue; end
        [az, el] = srp_block(xs{b}(i0:i1, :), M{b}, fs, cfg, N{b});   % front half-space only (planar board)
        ctr = mean(M{b}, 1);
        v = cfg.calib(3:5) - ctr;
        aztrue = atan2d(v(1), v(2));
        dz = wrap180(aztrue - az);
        fprintf('%-7s measured az %+6.1f  true %+6.1f  -> yaw error %+5.1f deg', cfg.names{b}, az, aztrue, dz);
        if abs(dz) > 30
            fprintf('   NOT corrected (> 30 deg: check facing / roll / which side the mics face)');
        elseif cfg.autoyaw
            Rz = [cosd(-dz) -sind(-dz) 0; sind(-dz) cosd(-dz) 0; 0 0 1];   % +dz: rotate toward +X
            M{b} = ((Rz * (M{b} - ctr)')' + ctr);
            N{b} = (Rz * N{b}')';
            fprintf('   (corrected)');
        end
        fprintf('\n');
    end
end

% ---------------------------------------------------------------- localisation
fprintf('\n=== 5. SRP-PHAT, band %d-%d Hz, block %.2f s ===\n', cfg.band, cfg.block);
i0 = max(1, round(cfg.seg(1) * fs) + 1); i1 = min(LA, round(cfg.seg(2) * fs));
blk = round(cfg.block * fs);
starts = i0:blk:(i1 - blk + 1);
sets = {1, 1:nb};
setname = {cfg.names{1}, 'ALL combined'};
res.t = (starts - 1 + blk / 2) / fs;
for s = 1:2
    Xs = [xs{sets{s}}]; Ms = vertcat(M{sets{s}});
    az = zeros(size(starts)); el = az;
    for k = 1:numel(starts)
        if s == 1, nrm = N{1}; else, nrm = []; end          % one planar board: front half only
        [az(k), el(k)] = srp_block(Xs(starts(k):starts(k) + blk - 1, :), Ms, fs, cfg, nrm);
    end
    res.az{s} = az; res.el{s} = el; res.set{s} = setname{s};
    fprintf('%-13s (%2d mics, aperture %.2f m): %d blocks\n', setname{s}, size(Ms, 1), aperture(Ms), numel(az));
end

% ---------------------------------------------------------------- errors vs truth
if isfield(cfg, 'truth') && ~isempty(cfg.truth)
    fprintf('\n=== 6. error vs known directions ===\n');
    fprintf('%-13s %-14s %8s %8s %10s\n', 'set', 'window (s)', 'med az', 'med el', 'within 10');
    for s = 1:2
        for w = 1:size(cfg.truth, 1)
            in = res.t >= cfg.truth(w, 1) & res.t <= cfg.truth(w, 2);
            if ~any(in), continue; end
            daz = abs(wrap180(res.az{s}(in) - cfg.truth(w, 3)));
            delv = abs(res.el{s}(in) - cfg.truth(w, 4));
            ok = mean(sqrt(daz .^ 2 + delv .^ 2) <= 10) * 100;
            fprintf('%-13s %5.1f-%-7.1f %7.1f° %7.1f° %8.0f %%\n', res.set{s}, cfg.truth(w, 1:2), median(daz), median(delv), ok);
        end
    end
end

% ---------------------------------------------------------------- plots
if cfg.plot
    figure('Name', 'tri_tdoa_combine', 'Position', [60 60 1200 700]);
    subplot(2, 2, 1); hold on; grid on; axis equal;
    cols = {[0.12 0.47 0.71], [1 0.5 0.05], [0.58 0.4 0.74]};
    for b = 1:nb, plot(M{b}(:, 1), M{b}(:, 2), 'o', 'Color', cols{b}, 'MarkerFaceColor', cols{b}); end
    plot(cfg.spk(1), cfg.spk(2), 'kp', 'MarkerSize', 12);
    legend([cfg.names, {'sync speaker'}], 'Location', 'best'); xlabel('x (m)'); ylabel('y (m)');
    title('mic positions (top view)');
    subplot(2, 2, 2); hold on; grid on;
    for s = 1:2, plot(res.t, res.az{s}, '.-'); end
    ylabel('azimuth (deg)'); xlabel('time (s)'); legend(res.set, 'Location', 'best'); title('azimuth');
    subplot(2, 2, 4); hold on; grid on;
    for s = 1:2, plot(res.t, res.el{s}, '.-'); end
    ylabel('elevation (deg)'); xlabel('time (s)'); title('elevation');
    subplot(2, 2, 3); hold on; grid on;
    for b = 2:nb
        [pa, pb] = pair_chirps(e{1}, e{b}, cfg.chirp.period * fs);
        p = polyfit(pa, pb - pa, 1);
        plot(pa / fs, ((pb - pa) - polyval(p, pa)) / fs * 1e6, 'o-');
    end
    ylabel('sync residual (us)'); xlabel('time (s)'); legend(cfg.names(2:end)); title('sync quality (|res| < 20 us = OK)');
end
res.M = M;
end

% =====================================================================================
function loc = odas_geometry(file, ch)
% mic coordinates from an ODAS cfg ('mu = (x, y, z)' in order), rows selected by channel
txt = fileread(file);
tok = regexp(txt, 'mu\s*=\s*\(\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*\)', 'tokens');
P = cellfun(@(t) str2double(t), tok, 'UniformOutput', false);
P = vertcat(P{:});
if isempty(P), error('%s: no "mu = (x, y, z)" mic positions found', file); end
mp = regexp(txt, 'map\s*:\s*\(([^)]*)\)', 'tokens', 'once');
if ~isempty(mp)
    map = str2num(mp{1}); %#ok<ST2NM>
else
    map = 1:size(P, 1);
end
loc = zeros(numel(ch), 3);
for i = 1:numel(ch)
    r = find(map == ch(i), 1);
    if isempty(r), error('%s: channel %d not in the mic list', file, ch(i)); end
    loc(i, :) = P(r, :);
end
end

function R = board_rotation(facing, roll, mount)
% board frame (x, y in board plane, z = normal) -> world (X right, Y forward, Z up)
n = [sind(facing), cosd(facing), 0];
if strcmp(mount, 'flat')
    r = [cosd(facing), -sind(facing), 0]; up = cross([0 0 1], r);
    R = [r', up', [0; 0; 1]];
else
    up = [0 0 1]; r = cross(up, n);                    % proper rotation (det +1)
    R = [r', up', n'];
end
Rr = [cosd(roll) -sind(roll) 0; sind(roll) cosd(roll) 0; 0 0 1];   % in-plane rotation of the board
R = R * Rr;
end

function arr = detect_chirps(x, tpl, fs, period)
% matched filter -> one arrival (sub-sample) per chirp
N = 2 ^ nextpow2(numel(x) + numel(tpl));
r = real(ifft(fft(x, N) .* conj(fft(tpl, N))));
r = r(1:numel(x));
X = fft(r); h = zeros(size(X)); h(1) = 1; h(2:floor(N / 2)) = 2; h(floor(N / 2) + 1) = 1;
env = abs(ifft(X .* h));                                % analytic envelope
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
        [~, q] = max(r(lo:hi)); q = q + lo - 1;          % carrier peak near the envelope max
        a = r(q - 1); b = r(q); c = r(q + 1);
        arr(end + 1, 1) = (q - 1) + 0.5 * (a - c) / (a - 2 * b + c); %#ok<AGROW>
        i = j + 1;
    else
        i = i + 1;
    end
end
end

function [pa, pb] = pair_chirps(ea, eb, per)
% pair chirps of board A and board B (offset between recordings < period / 2)
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

function [az, el] = srp_block(X, Mpos, fs, cfg, nrm)
% SRP-PHAT over a direction grid (coarse, then fine around the peak); far field unless cfg.range
% nrm: board normal -> search only in front of that board (a planar array cannot tell front/back)
nf = cfg.nfft; hop = nf / 2;
w = 0.5 - 0.5 * cos(2 * pi * (0:nf - 1)' / nf);
f = (0:nf / 2)' * fs / nf;
sel = find(f >= cfg.band(1) & f <= cfg.band(2));
fr = f(sel);
nfr = floor((size(X, 1) - nf) / hop) + 1;
Y = zeros(numel(sel), size(X, 2), max(nfr, 1));
for k = 1:nfr
    F = fft(X((k - 1) * hop + (1:nf), :) .* w);
    F = F(sel, :);
    Y(:, :, k) = F ./ (abs(F) + 1e-12);                % PHAT
end
[AZ, EL] = meshgrid(cfg.grid_az, cfg.grid_el);
AZ = AZ(:); EL = EL(:);
if isfield(cfg, 'upper_only') && cfg.upper_only     % sources above the table plane
    AZ = AZ(EL >= 0); EL = EL(EL >= 0);
end
if ~isempty(nrm)
    keep = [sind(AZ) .* cosd(EL), cosd(AZ) .* cosd(EL), sind(EL)] * nrm(:) > 0.05;
    AZ = AZ(keep); EL = EL(keep);
end
P = srp_eval(Y, Mpos, fr, AZ(:), EL(:), cfg);
[~, i] = max(P);
az0 = AZ(i); el0 = EL(i);
[AZ, EL] = meshgrid(az0 + (-4:0.5:4), max(-90, min(90, el0 + (-6:0.5:6))));
P = srp_eval(Y, Mpos, fr, AZ(:), EL(:), cfg);
[~, i] = max(P);
az = wrap180(AZ(i)); el = EL(i);
end

function P = srp_eval(Y, Mpos, fr, az, el, cfg)
c0 = 343;
u = [sind(az) .* cosd(el), cosd(az) .* cosd(el), sind(el)];         % D x 3
if isfield(cfg, 'range') && isfinite(cfg.range)
    ctr = mean(Mpos, 1);
    S = ctr + cfg.range * u;                                         % source points
    tau = zeros(size(u, 1), size(Mpos, 1));
    for m = 1:size(Mpos, 1)
        tau(:, m) = sqrt(sum((S - Mpos(m, :)) .^ 2, 2)) / c0;
    end
else
    tau = -(u * Mpos') / c0;                                         % D x M arrival delays
end
P = zeros(size(u, 1), 1);
for fi = 1:numel(fr)
    A = exp(1j * 2 * pi * fr(fi) * tau);                             % undo each mic delay
    Yf = squeeze(Y(fi, :, :));                                       % M x frames
    if size(Y, 2) == 1, Yf = Yf(:)'; end
    P = P + sum(abs(A * Yf) .^ 2, 2);
end
end

function a = wrap180(a)
a = mod(a + 180, 360) - 180;
end

function d = aperture(M)
d = 0;
for i = 1:size(M, 1)
    d = max(d, max(sqrt(sum((M - M(i, :)) .^ 2, 2))));
end
end

function [R, t, rms_mm] = kabsch(L, W)
% best rotation R + translation t with W ~ R*L + t  (rows = points)
cl = mean(L, 1); cw = mean(W, 1);
H = (L - cl)' * (W - cw);
[U, ~, V] = svd(H);
D = diag([1 1 sign(det(V * U'))]);
R = V * D * U';
t = cw - (R * cl')';
E = (R * L')' + t - W;
rms_mm = 1000 * sqrt(mean(sum(E .^ 2, 2)));
end

function r = ternary(c, a, b)
if c, r = a; else, r = b; end
end
