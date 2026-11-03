function tri_level_check(seconds)
% TRI_LEVEL_CHECK  Live sound meter of the 3 boards: is the phone heard?
%   >> tri_level_check        % 60 s (close the figure to stop earlier)
% Green bar (above -60 dBFS) = loud enough for calibration / localisation.
if nargin < 1, seconds = 60; end
fs = 48000; frame = round(0.1 * fs);
dev = tri_devices();
rd = cell(1, 3);
for b = 1:3
    rd{b} = audioDeviceReader('Driver', dev(b).driver, 'Device', dev(b).name, ...
        'NumChannels', dev(b).nch, 'SampleRate', fs, 'SamplesPerFrame', frame, 'OutputDataType', 'double');
    setup(rd{b});
end
cleanup = onCleanup(@() cellfun(@release, rd)); %#ok<NASGU>
names = {'UMA16A', 'UMA8', 'UMA16B'};
fig = figure('Name', 'tri_level_check - close to stop', 'Position', [200 200 600 420]);
ax = axes(fig); hold(ax, 'on'); grid(ax, 'on');
h = bar(ax, 1:3, [-100 -100 -100], 'BaseValue', -100, 'FaceColor', 'flat');
yline(ax, -60, 'g--', 'loud enough', 'LineWidth', 1.5);
yline(ax, -65, 'r--', 'cfg.live\_gate');
set(ax, 'XTick', 1:3, 'XTickLabel', names); ylim(ax, [-100 0]); ylabel(ax, 'level (dBFS)');
for k = 1:10, for b = 1:3, rd{b}(); end, end
t0 = tic;
while ishghandle(fig) && toc(t0) < seconds
    lev = zeros(1, 3);
    for b = 1:3
        x = rd{b}();
        x = x(:, any(x, 1));                          % live channels only
        lev(b) = 20 * log10(rms(x(:)) + eps);
    end
    h.YData = lev;
    h.CData = repmat([0.85 0.2 0.2], 3, 1);
    h.CData(lev > -60, :) = repmat([0.2 0.7 0.2], sum(lev > -60), 1);
    title(ax, sprintf('UMA16A %.0f   UMA8 %.0f   UMA16B %.0f dBFS', lev));
    drawnow limitrate;
end
end
