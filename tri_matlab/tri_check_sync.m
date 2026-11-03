function tri_check_sync(seconds)
% TRI_CHECK_SYNC  Quick test: record a few seconds with the chirp playing and check
% that every array hears one chirp per second.  Run before a real tri_record.
%
%   >> tri_check_sync          % 10 s
if nargin < 1, seconds = 10; end
tri_record(seconds);
fs = 48000;
tpl = sync_chirp_template(fs, 6000, 9000, 0.03);
files = {'uma16a.wav', 'uma8.wav', 'uma16b.wav'};
ok = true;
fprintf('\n=== sync chirp check (expect ~%d chirps, 1 per second) ===\n', seconds);
for b = 1:3
    x = audioread(files{b});
    x = x(:, 1) - mean(x(:, 1));
    r = abs(hilbert(filter(flipud(tpl), 1, x)));        % matched-filter envelope
    nsec = floor(numel(r) / fs);
    pk = zeros(1, nsec); snr = zeros(1, nsec);
    for s = 1:nsec
        seg = r((s - 1) * fs + 1:s * fs);
        [m, i] = max(seg);
        pk(s) = mod(i, fs); snr(s) = 20 * log10(m / median(seg));
    end
    good = snr > 20;                                    % clear chirp in that second
    jitter = std(pk(good)) / fs * 1e3;                  % chirps must repeat at the same phase
    pass = sum(good) >= 0.8 * nsec && jitter < 1;
    ok = ok && pass;
    fprintf('%-11s %2d/%2d seconds with a clear chirp, median SNR %4.1f dB, timing jitter %.2f ms  -> %s\n', ...
        files{b}, sum(good), nsec, median(snr), jitter, ternary(pass, 'PASS', 'FAIL'));
end
if ok
    fprintf('\nSync OK. Now: tri_record(60) while making your test sound, then tri_config\n');
else
    fprintf(['\nFAIL: the arrays do not hear the chirp. Unmute the laptop speaker, volume ~80-100 %%,\n' ...
             'unplug headphones, and keep the laptop where cfg.spk says.\n']);
end
end

function r = ternary(c, a, b)
if c, r = a; else, r = b; end
end
