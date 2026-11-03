function [det, spec] = freq_detect(X, fs, F, blacklist)
% FREQ_DETECT  Find a drone's harmonic comb (fundamental f0 + harmonics), allowing f0 to drift.
%
%   X          1 x nb cell, samples x mics of every board (the last cfg.freq.listen_sec seconds)
%   F          cfg.freq (see tri_config.m)
%   blacklist  frequencies (Hz) of known stationary electrical tones, ignored (optional)
%   det   .f .snr_db     bin frequencies (0..fmax) and the averaged whitened spectrum (dB above noise)
%         .detected      true if a comb passes every test;  .reason  why not (text)
%         .f0            mean fundamental (Hz);  .f0_track  f0 in every frame (Hz)
%         .harm_h        orders of the harmonics that count;  .harm_f  their mean frequency (Hz)
%         .harm_bins     their bin in every frame (harmonics x frames);  .harm_snr  SNR (dB)
%         .harm_coh      drift-compensated steadiness (0..1);  .score / .thr  comb score / threshold
%   spec  1 x nb cell, bins x mics x frames complex spectra (reused for the direction)
%
% Method:
%   1. Welch frames (nfft, Hann, 50 % overlap), power averaged over all 25 mics, whitened per frame
%      (SNR = power / running median over +-100 Hz).
%   2. Per frame, every f0 of a log-spaced grid gets the sum of its harmonics' SNR (h*f0 <= fmax).
%   3. Dynamic programming over the frames: the f0 path with the highest total score, f0 changing by
%      at most max_rate per frame (motor speed changes).  Harmonic h follows h*f0(t), i.e. it drifts
%      in proportion to its order.
%   4. The best paths are checked in order; the first that passes all tests is the detection:
%      - harmonic h counts if, along its own track, it is peak_db above the noise, min_contrast above
%        the bins 18-24 Hz away (narrow: rejects room / box / speaker resonances), and steady:
%        the frame-to-frame phase ADVANCE changes smoothly (coherence of its change >= min_coh_h).
%        This accepts a drifting tone but not noise passed through a resonance, whose phase is random;
%      - the fundamental must carry its own structure: at least 2 of its harmonics 1-3 count
%        (a sub-multiple such as f0/11 has none of them and is rejected);
%      - at least min_harm harmonics count, their mean steadiness >= min_coh, score >= min_score.
nf = F.nfft; hop = nf / 2;
w = 0.5 - 0.5 * cos(2 * pi * (0:nf - 1)' / nf);
f = (0:nf / 2)' * fs / nf;
keep = f <= F.fmax;
f = f(keep); nbin = numel(f); df = f(2) - f(1);
nb = numel(X);
spec = cell(1, nb); Pk = 0; nm = 0;
for b = 1:nb
    x = X{b} - mean(X{b}, 1);
    nfr = floor((size(x, 1) - nf) / hop) + 1;
    S = zeros(nbin, size(x, 2), nfr);
    for k = 1:nfr
        Z = fft(x((k - 1) * hop + (1:nf), :) .* w);
        S(:, :, k) = Z(keep, :);
    end
    spec{b} = S;
    Pk = Pk + squeeze(sum(abs(S) .^ 2, 2));
    nm = nm + size(x, 2);
end
Pk = Pk / nm; K = size(Pk, 2);
nw = 2 * round(100 / df) + 1;
noise = movmedian(Pk, nw, 1);
snrk = 10 * log10(Pk ./ max(noise, eps));             % bins x frames
P = mean(Pk, 2);
snr_db = 10 * log10(P ./ max(movmedian(P, nw), eps));
bl = false(nbin, 1);
if nargin >= 4
    for fb = blacklist(:)', bl = bl | abs(f - fb) <= 1.5 * df; end
end
det = struct('f', f, 'snr_db', snr_db, 'detected', false, 'reason', '', 'f0', NaN, 'f0_track', [], ...
    'harm_h', [], 'harm_f', [], 'harm_bins', [], 'harm_snr', [], 'harm_coh', [], 'score', NaN, 'thr', NaN);
if K < 4, det.reason = 'too few frames'; return; end

% 2. per-frame harmonic score on a log-spaced f0 grid (relative step f0_step_rel)
f0s = exp(log(F.f0_range(1)):log(1 + F.f0_step_rel):log(F.f0_range(2)))';
h = 1:F.n_harm;
hf = f0s * h; valid = hf <= F.fmax - df;
idx = min(round(hf / df) + 1, nbin);
c = min(max(movmax(snrk, 3, 1) - 3, 0), 17);          % dB above 3 dB, best bin within +-1 bin
c(bl, :) = 0;
H = zeros(numel(f0s), K);
for k = 1:K
    ck = c(:, k); v = ck(idx); v(~valid) = 0;
    H(:, k) = sum(v, 2);
end
% 3. dynamic programming: f0 may change by max_rate per frame
wn = max(1, round(log(1 + F.max_rate) / log(1 + F.f0_step_rel)));
D = H(:, 1); arg = zeros(numel(f0s), K, 'uint32');
for k = 2:K
    [m, a] = winmax(D, wn);
    D = H(:, k) + m; arg(:, k) = a;
end
det.thr = max(median(D) + F.cfar_k * 1.4826 * mad(D, 1), F.min_score);
cand = find(islocalmax(D) & D > det.thr);
[~, o] = sort(D(cand), 'descend');
cand = cand(o(1:min(end, F.n_cand)));
if isempty(cand), det.reason = 'no comb above the threshold'; det.score = max(D); return; end

% 4. check the best paths in order
det.reason = 'no candidate passed the harmonic tests';
for i = cand'
    p = zeros(1, K); p(K) = i;
    for k = K:-1:2, p(k - 1) = arg(p(k), k); end
    f0k = f0s(p)';
    r = check_path(f0k);
    if r.ok
        det.detected = true; det.reason = '';
        det.f0 = mean(f0k); det.f0_track = f0k;
        det.harm_h = r.h; det.harm_f = r.f; det.harm_bins = r.bins;
        det.harm_snr = r.snr; det.harm_coh = r.coh; det.score = r.score;
        return
    end
end

    function r = check_path(f0k)
        hs = find(F.n_harm >= h & h * max(f0k) <= F.fmax - 5 * df);
        n = numel(hs);
        bins = zeros(n, K); snrh = zeros(1, n); con = zeros(1, n); coh = zeros(1, n);
        for q = 1:n
            bq = round(hs(q) * f0k / df) + 1;
            for k2 = 1:K                               % strongest bin within +-1 bin, per frame
                lo = max(1, bq(k2) - 1); [~, j] = max(Pk(lo:min(nbin, bq(k2) + 1), k2)); bq(k2) = lo + j - 1;
            end
            bins(q, :) = bq;
            lin = sub2ind(size(Pk), bq, 1:K);
            snrh(q) = 10 * log10(mean(Pk(lin)) / mean(noise(lin)));
            sidev = zeros(1, K);
            for k2 = 1:K
                sb = bq(k2) + [-4 -3 3 4]; sb = sb(sb >= 1 & sb <= nbin);
                sidev(k2) = max(Pk(sb, k2));
            end
            con(q) = 10 * log10(mean(Pk(lin)) / max(mean(sidev), eps));
            % steadiness: Z_k = sum over all mics of X_k X*_(k-2) (phase advance over 2 frames);
            % a tone (even drifting) changes this advance smoothly, resonance noise randomly
            Zk = zeros(1, K - 2);
            for k2 = 3:K
                for bb = 1:nb
                    Zk(k2 - 2) = Zk(k2 - 2) + sum(spec{bb}(bq(k2), :, k2) .* conj(spec{bb}(bq(k2 - 2), :, k2 - 2)));
                end
            end
            dz = Zk(2:end) .* conj(Zk(1:end - 1));
            coh(q) = abs(sum(dz)) / max(sum(abs(dz)), eps);
        end
        inbl = false(1, n);
        for q = 1:n, inbl(q) = any(bl(bins(q, :))); end
        counts = snrh >= F.peak_db & con >= F.min_contrast & coh >= F.min_coh_h & ~inbl;
        low = find(hs <= 3);
        r.ok = numel(low) >= 1 && sum(counts(low)) >= min(2, numel(low)) ...   % own structure
            && sum(counts) >= F.min_harm && mean(coh(counts)) >= F.min_coh;
        r.score = sum(min(snrh(counts), 20));
        r.ok = r.ok && r.score >= F.min_score;
        r.h = hs(counts); r.f = mean(f(bins(counts, :)), 2)';
        r.bins = bins(counts, :); r.snr = snrh(counts); r.coh = coh(counts);
    end
end

% =====================================================================================
function [m, a] = winmax(D, w)
% running maximum of D over [i-w, i+w] and the index where it is reached
n = numel(D); m = D; a = uint32(1:n)';
for s = [-w:-1, 1:w]
    j = (1:n)' + s; ok = j >= 1 & j <= n;
    v = -inf(n, 1); v(ok) = D(j(ok));
    better = v > m; m(better) = v(better); a(better) = j(better);
end
end
