function sync_chirp_make(outfile, minutes)
% SYNC_CHIRP_MAKE  Sync track for the 3-array recording.
%   One 30 ms linear chirp 6-9 kHz every 1.0 s (above the 200-4000 Hz analysis band,
%   so it never disturbs the drone/phone signal).  Play it from the SYNC SPEAKER
%   (laptop speaker) during the whole recording.
%
%   sync_chirp_make('sync_chirp.wav', 6)
if nargin < 1, outfile = 'sync_chirp.wav'; end
if nargin < 2, minutes = 6; end
fs = 48000; dur = 0.03; f0 = 6000; f1 = 9000; period = 1.0;
c = sync_chirp_template(fs, f0, f1, dur);
n = round(minutes * 60 / period);
y = zeros(round(n * period * fs), 1);
for k = 0:n-1
    i0 = round(k * period * fs) + 1;
    y(i0:i0 + numel(c) - 1) = c;
end
audiowrite(outfile, 0.8 * y, fs);
fprintf('wrote %s: %d chirps, %.0f s\n', outfile, n, numel(y) / fs);
end
