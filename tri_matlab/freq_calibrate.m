function freq_calibrate(cfg)
% FREQ_CALIBRATE  Board-orientation calibration for the frequency mode (same method as tdoa_live).
%
%   The phone plays WHITE NOISE at the points cfg.tdoa.cal (move_sec to move, rec_sec still).
%   Board orientations come from tdoa_calibrate (mounting priors cfg.tdoa.prior), the systematic
%   error from tdoa_systematic, exactly as in tdoa_live.  Saved to tdoa_calibration.mat together
%   with the board positions used (cfg.tri.pos), so a later run can tell when they changed.
%   Called by freq_live when cfg.freq.calibrate = true.
fs = cfg.fs; frame = round(0.1 * fs); blk = round(cfg.block * fs);
S = cfg.tdoa;
B = tdoa_boards(cfg);
cal = S.cal; if size(cal, 2) == 2, cal(:, 3) = S.src_z; end
G = 1; if isfield(cfg, 'mic_gain'), G = cfg.mic_gain; end
here = fileparts(mfilename('fullpath'));

fig = figure('Name', 'Calibration (frequency mode)', 'Color', 'w', 'Position', [200 200 700 520]);
ax = axes(fig); hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
P = vertcat(cfg.tri.pos{:});
plot(ax, P(:, 1), P(:, 2), 'ks', 'MarkerSize', 10, 'MarkerFaceColor', [0.5 0.5 0.5]);
plot(ax, cal(:, 1), cal(:, 2), 'kx', 'MarkerSize', 12, 'LineWidth', 2);
text(ax, cal(:, 1) + 0.04, cal(:, 2) + 0.06, compose('cal %d (z=%.2f m)', [1:size(cal, 1); cal(:, 3)']'));
hNow = plot(ax, nan, nan, 'ro', 'MarkerSize', 26, 'LineWidth', 3);
lim = max(abs([cal(:, 1); cal(:, 2)])) + 0.4; xlim(ax, [-lim lim]); ylim(ax, [-lim lim]);
xlabel(ax, 'X (m)  right'); ylabel(ax, 'Y (m)  toward UMA-8');

dev = tri_devices();
rd = cell(1, 3);
for b = 1:3
    rd{b} = audioDeviceReader('Driver', dev(b).driver, 'Device', dev(b).name, ...
        'NumChannels', dev(b).nch, 'SampleRate', fs, 'SamplesPerFrame', frame, 'OutputDataType', 'double');
    setup(rd{b});
end
cleanup = onCleanup(@() cellfun(@release, rd)); %#ok<NASGU>
for k = 1:10, for b = 1:3, rd{b}(); end, end         % warm-up

taucal = cell(size(cal, 1), 3);
for p = 1:size(cal, 1)
    set(hNow, 'XData', cal(p, 1), 'YData', cal(p, 2));
    T = cell(1, 3);
    for phase = 1:2
        dur = [S.move_sec, S.rec_sec]; dur = dur(phase);
        buf = {[], [], []}; t0 = tic;
        while toc(t0) < dur && ishghandle(fig)
            for b = 1:3
                [x, lost] = rd{b}();
                buf{b} = [buf{b}; zeros(double(lost), numel(cfg.use_ch{b})); G * x(:, cfg.use_ch{b})];
            end
            if phase == 1
                title(ax, sprintf('MOVE the phone (white noise) to cal %d - %.0f s', p, dur - toc(t0)), 'FontSize', 13);
            else
                title(ax, sprintf('cal %d: keep the phone STILL - %.0f s', p, dur - toc(t0)), 'FontSize', 13);
            end
            drawnow limitrate;
            if phase == 2 && all(cellfun(@(z) size(z, 1), buf) >= blk)
                Xb = cellfun(@(z) z(end - blk + 1:end, :), buf, 'UniformOutput', false);
                buf = {[], [], []};
                if toc(t0) < S.settle_sec || 20 * log10(rms(Xb{1}(:)) + eps) < cfg.live_gate, continue; end
                for b = 1:3
                    T{b}(:, end + 1) = tdoa_estimate(Xb{b} - mean(Xb{b}, 1), B(b).pairs, B(b).maxlag, fs, S);
                end
            elseif phase == 1
                buf = {[], [], []};
            end
        end
        if ~ishghandle(fig), error('Calibration cancelled (figure closed).'); end
    end
    if isempty(T{1})
        error('cal %d: the phone was not heard (nothing louder than cfg.live_gate = %d dBFS).', p, cfg.live_gate);
    end
    for b = 1:3, taucal{p, b} = median(T{b}, 2); end
    fprintf('cal %d done (%d blocks)\n', p, size(T{1}, 2));
end
for b = 1:3
    Tb = [taucal{:, b}];
    [B(b).pose, e] = tdoa_calibrate(B(b).Mloc, B(b).centre, B(b).pairs, Tb, cal, S.src_z, S.flat(b), prior_of(S, b));
    tdoa_print_pose(B(b), e);
end
sys = tdoa_systematic(B, taucal, cal, cfg);
pose = [B.pose]; when = char(datetime('now')); pos = cfg.tri.pos; %#ok<NASGU>
save(fullfile(here, 'tdoa_calibration.mat'), 'pose', 'when', 'sys', 'pos');
fprintf('Calibration saved (with the board positions) to %s\n', fullfile(here, 'tdoa_calibration.mat'));
close(fig);
end
