function tdoa_session(cfg)
% TDOA_SESSION  Guided recording + evaluation of the own TDOA method.
%
%   The figure tells you where to put the phone (white noise, full volume):
%     cal 1..n   cfg.tdoa.cal    -> 3-D calibration of every board
%     test 1..m  cfg.tdoa.test   -> known positions used to MEASURE the error
%   Everything is recorded (live_sessions/<stamp>_*.wav), then tdoa_offline runs and
%   the error at every test point is printed and saved (<stamp>_tdoa_eval.csv).
%   Called from tri_config.m when cfg.mode = 'tdoa_session'.
if nargin < 1                                        % typed alone: take the settings from tri_config
    fprintf('No settings given - running tri_config (set cfg.mode = ''tdoa_session'' there).\n');
    tri_config;
    return
end
fs = cfg.fs; frame = round(0.1 * fs);
S = cfg.tdoa;
steps = [S.cal; S.test];
nCal = size(S.cal, 1);
label = [compose('cal %d', 1:nCal), compose('test %d', 1:size(S.test, 1))];

fig = figure('Name', 'tdoa_session', 'Position', [100 100 900 600]);
ax = axes(fig); hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
P = vertcat(cfg.tri.pos{:});
plot(ax, P(:, 1), P(:, 2), 'ks', 'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.5 0.5]);
text(ax, P(:, 1) + 0.03, P(:, 2) - 0.05, cfg.names);
plot(ax, steps(:, 1), steps(:, 2), 'bx', 'MarkerSize', 10, 'LineWidth', 2);
text(ax, steps(:, 1) + 0.03, steps(:, 2) + 0.04, label);
hNow = plot(ax, nan, nan, 'ro', 'MarkerSize', 22, 'LineWidth', 3);
lim = max(abs([steps(:, 1); steps(:, 2)])) + 0.4;
xlim(ax, [-lim lim]); ylim(ax, [-lim lim]);
xlabel(ax, 'X (m)  right'); ylabel(ax, 'Y (m)  toward UMA-8');
drawnow;

dev = tri_devices();
outdir = fullfile(fileparts(mfilename('fullpath')), 'live_sessions');
if ~isfolder(outdir), mkdir(outdir); end
stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
rd = cell(1, 3); wr = cell(1, 3);
for b = 1:3
    rd{b} = audioDeviceReader('Driver', dev(b).driver, 'Device', dev(b).name, ...
        'NumChannels', dev(b).nch, 'SampleRate', fs, 'SamplesPerFrame', frame, 'OutputDataType', 'double');
    setup(rd{b});
    wr{b} = dsp.AudioFileWriter(fullfile(outdir, [stamp '_' dev(b).file]), 'SampleRate', fs, 'DataType', 'int32');
end
cleanup = onCleanup(@() cellfun(@release, [wr rd])); %#ok<NASGU>
for k = 1:10, for b = 1:3, rd{b}(); end, end         % warm-up

win = zeros(size(steps, 1), 2);
nRead = 0;                                            % frames written = recording time base
for s = 1:size(steps, 1)
    set(hNow, 'XData', steps(s, 1), 'YData', steps(s, 2));
    for phase = 1:2
        dur = [S.move_sec, S.rec_sec]; dur = dur(phase);
        n0 = nRead;
        while (nRead - n0) * frame / fs < dur
            for b = 1:3
                [x, lost] = rd{b}();
                tri_write(wr{b}, x, lost);
            end
            nRead = nRead + 1;
            left = dur - (nRead - n0) * frame / fs;
            if phase == 1
                msg = sprintf('MOVE the phone to %s  (x=%.2f, y=%.2f, z=%.2f m), white noise ON - %.0f s', label{s}, steps(s, :), left);
            else
                msg = sprintf('%s: keep the phone STILL - recording %.0f s', upper(label{s}), left);
            end
            if mod(nRead, 5) == 0, title(ax, msg, 'FontSize', 14); drawnow limitrate; end
        end
        if phase == 2
            win(s, :) = [n0 * frame / fs + S.settle_sec, nRead * frame / fs];
            fprintf('%-7s recorded (%.1f-%.1f s)\n', label{s}, win(s, :));
        end
    end
end
clear cleanup                                        % close the files before analysing
title(ax, 'Recording done - analysing ...', 'FontSize', 14); drawnow;

tdoa_evaluate(cfg, stamp, win, ax);
end
