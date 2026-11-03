function c = sync_chirp_template(fs, f0, f1, dur)
% Linear chirp f0 -> f1 Hz, dur seconds, Hann-tapered (used to make AND to detect the sync chirps)
t = (0:round(dur * fs) - 1)' / fs;
c = sin(2 * pi * (f0 * t + (f1 - f0) / (2 * dur) * t .^ 2));
w = 0.5 - 0.5 * cos(2 * pi * (0:numel(t) - 1)' / (numel(t) - 1));
c = c .* w;
end
