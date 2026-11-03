function res = tdoa_offline(cfg, stamp, calwin)
% TDOA_OFFLINE  Run the TDOA method on a saved live session (live_sessions/<stamp>_*.wav).
%
%   calwin  n x 2 [t0 t1] seconds when the phone was at the calibration points cfg.tdoa.cal
%   res     .t .pos [x y z] .std .rms_us  per block, .pose  per board
%   Writes live_sessions/<stamp>_tdoa_offline.csv
fs = cfg.fs;
dir0 = fullfile(fileparts(mfilename('fullpath')), 'live_sessions');
B = tdoa_boards(cfg);
files = {'uma16a.wav', 'uma8.wav', 'uma16b.wav'};
x = cell(1, 3);
for b = 1:3
    y = audioread(fullfile(dir0, [stamp '_' files{b}]), 'native');   % int32: half the memory
    x{b} = double(y(:, B(b).ch)) / 2 ^ 31;
    if isfield(cfg, 'mic_gain'), x{b} = cfg.mic_gain * x{b}; end   % same amplifier as live
end
blk = round(cfg.block * fs);
n = min(cellfun(@(z) size(z, 1), x));
starts = 1:blk:(n - blk + 1);
t = (starts - 1 + blk / 2) / fs;

% TDOAs of every block, every board
tau = cell(numel(starts), 3); conf = tau; lev = zeros(numel(starts), 1);
for k = 1:numel(starts)
    for b = 1:3
        X = x{b}(starts(k):starts(k) + blk - 1, :);
        X = X - mean(X, 1);
        [tau{k, b}, conf{k, b}] = tdoa_estimate(X, B(b).pairs, B(b).maxlag, fs, cfg.tdoa);
    end
    lev(k) = 20 * log10(rms(reshape(x{1}(starts(k):starts(k) + blk - 1, :), [], 1)) + eps);
end

% calibration: median TDOAs inside each window -> 3-D pose of each board
ncal = size(calwin, 1);
for b = 1:3
    taucal = zeros(size(B(b).pairs, 1), ncal);
    for p = 1:ncal
        in = t >= calwin(p, 1) & t <= calwin(p, 2) & lev' >= cfg.live_gate;
        if ~any(in)
            error(['cal %d (%.0f-%.0f s): no block louder than cfg.live_gate = %d dBFS (level there: %.0f dBFS).\n' ...
                   'The phone was not heard: white noise at FULL volume, speaker facing the boards.'], ...
                p, calwin(p, :), cfg.live_gate, max(lev(t >= calwin(p, 1) & t <= calwin(p, 2))));
        end
        taucal(:, p) = median([tau{in, b}], 2);
    end
    [B(b).pose, e] = tdoa_calibrate(B(b).Mloc, B(b).centre, B(b).pairs, taucal, cfg.tdoa.cal(1:ncal, :), cfg.tdoa.src_z, cfg.tdoa.flat(b), prior_of(cfg.tdoa, b));
    tdoa_print_pose(B(b), e);
end

% localisation of every block
res.t = t; res.pos = nan(numel(starts), 3); res.std = res.pos; res.rms_us = nan(numel(starts), 1);
res.level = lev; res.cov = nan(numel(starts), 3, 3);
for k = 1:numel(starts)
    if lev(k) < cfg.live_gate, continue; end
    [res.pos(k, :), info] = tdoa_locate(B, tau(k, :), conf(k, :), cfg);
    res.std(k, :) = info.std_xy; res.rms_us(k) = info.rms_us; res.cov(k, :, :) = info.cov;
end
res.pose = [B.pose];

out = fullfile(dir0, [stamp '_tdoa_offline.csv']);
fid = fopen(out, 'w');
fprintf(fid, 't_s,level_dBFS,x,y,z,std_x,std_y,std_z,rms_us\n');
fprintf(fid, '%.2f,%.1f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.1f\n', [t(:), lev, res.pos, res.std, res.rms_us]');
fclose(fid);
fprintf('wrote %s\n', out);
end
