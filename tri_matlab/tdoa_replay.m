function tdoa_replay(cfg, stamp, t_from)
% TDOA_REPLAY  Replay a saved session (live_sessions/<stamp>_*.wav) with the live display.
%
%   tdoa_replay(cfg, '20261102_092801')        % from the end of the calibration
%   tdoa_replay(cfg, '20261102_092801', 0)     % from the start of the recording
% Uses the calibration in tdoa_calibration.mat (the one made during that session if no
% calibration was made since).  Saves the final figure as live_sessions/<stamp>_view.png.
here = fileparts(mfilename('fullpath'));
dir0 = fullfile(here, 'live_sessions');
S = cfg.tdoa; fs = cfg.fs; blk = round(cfg.block * fs);
if nargin < 3, t_from = size(S.cal, 1) * (S.move_sec + S.rec_sec); end
B = tdoa_boards(cfg);
L = load(fullfile(here, 'tdoa_calibration.mat'));
for b = 1:3, B(b).pose = L.pose(b); end
sys = []; if isfield(L, 'sys'), sys = L.sys; end
fprintf('replaying %s with the calibration of %s\n', stamp, L.when);
files = {'uma16a.wav', 'uma8.wav', 'uma16b.wav'};
x = cell(1, 3);
for b = 1:3
    y = audioread(fullfile(dir0, [stamp '_' files{b}]), 'native');
    x{b} = double(y(:, B(b).ch)) / 2 ^ 31;
    if isfield(cfg, 'mic_gain'), x{b} = cfg.mic_gain * x{b}; end   % same amplifier as live
end
n = min(cellfun(@(z) size(z, 1), x));
V = tdoa_view('init', cfg, B);
tdoa_view('mics', V, B);
for s = round(t_from * fs) + 1:blk:n - blk + 1
    if ~ishghandle(V.fig), return; end
    t = (s - 1 + blk) / fs;
    X = cellfun(@(z) z(s:s + blk - 1, :), x, 'UniformOutput', false);
    lev = 20 * log10(rms(X{1}(:)) + eps);
    if lev < cfg.live_gate, continue; end
    tau = cell(1, 3); conf = tau;
    for b = 1:3
        [tau{b}, conf{b}] = tdoa_estimate(X{b} - mean(X{b}, 1), B(b).pairs, B(b).maxlag, fs, S);
    end
    [p, info] = tdoa_locate(B, tau, conf, cfg);
    C = tdoa_region(p, info.cov, sys);
    V = tdoa_view('update', V, t, p, C, info.rms_us < S.max_rms_us, lev);
end
png = fullfile(dir0, [stamp '_view.png']);
exportgraphics(V.fig, png, 'Resolution', 110);
fprintf('figure saved to %s\n', png);
end
