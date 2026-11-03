function [tau, conf, info] = tdoa_estimate(X, pairs, maxlag, fs, opt)
% TDOA_ESTIMATE  Time differences of arrival for the mic pairs of ONE board (one clock).
%
%   X       samples x mics, one block
%   pairs   P x 2 mic indices [i j];  tau(p) = arrival time at i - arrival time at j (s)
%   maxlag  P x 1 physical limit: |tau| <= mic spacing / c (s)
%   opt     .nfft   frame length for the cross-spectra
%           .band   [f0 f1] Hz used
%           .bands  K x 2 sub-bands (Hz) for the consensus
%           .up     lag-grid oversampling (lag step = 1 / (fs * up))
%           .noise  (optional) noise model of this board from tdoa_noise ([] = none)
%           .legacy (optional) true = exact previous behaviour (no noise model, equal bands)
%   conf    P x 1 confidence in [0, 1]
%   info    .f .Sxx .Sxy .nfr   spectra of this block (for tdoa_noise)
%           .R .lags            consensus correlation curve per pair (lag x pair), NOT masked
%                               (for tdoa_srp: all pairs vote before any decision)
%           .evid               1 x P evidence = sum over bins of the weight (SNR-like)
%           .snr_db             block power above the learnt noise floor (dB, NaN without model)
%
% Method (same idea as freq mode: only what stands above the noise counts):
%   1. Welch cross-spectra over the block.
%   2. The learnt NOISE cross-spectrum is subtracted: stationary sounds (fans, PC, electrical
%      tones that are coherent between mics) no longer pull the lag toward their own direction
%      or toward 0.
%   3. Signal coherence = |Sxy - Nxy|^2 / (Sxx Syy) minus the bias of a K-frame estimate
%      (1.1 / K): a noise-only bin gets weight ~0 instead of ~1/K.
%   4. Weight = PHAT x g2 / (1 - g2), evaluated on a fractional-lag grid inside the limit.
%   5. Sub-band consensus WEIGHTED by each band's evidence: a band with no signal no longer
%      gets the same vote as the band that carries the sound (it did before: each band was
%      normalised to 1, so at low SNR the noise-only bands added random peaks).
%   6. Parabolic refinement of the winning lag; confidence = consensus peak x mean coherence.
if ~isfield(opt, 'legacy') || isempty(opt.legacy), opt.legacy = false; end
noise = [];
if isfield(opt, 'noise'), noise = opt.noise; end

nf = opt.nfft; hop = nf / 2;
w = 0.5 - 0.5 * cos(2 * pi * (0:nf - 1)' / nf);
nfr = floor((size(X, 1) - nf) / hop) + 1;
if nfr < 2, error('tdoa_estimate: block shorter than 2 frames of %d samples', nf); end
f = (0:nf / 2)' * fs / nf;
sel = f >= opt.band(1) & f <= opt.band(2);
f = f(sel);

S = zeros(numel(f), size(X, 2), nfr);                 % spectra: freq x mic x frame
for k = 1:nfr
    F = fft(X((k - 1) * hop + (1:nf), :) .* w);
    S(:, :, k) = F(sel, :);
end
i = pairs(:, 1); j = pairs(:, 2);
P = size(pairs, 1);
Sxy = mean(S(:, i, :) .* conj(S(:, j, :)), 3);        % cross-spectra, freq x pair
Sxx = mean(abs(S) .^ 2, 3);                           % auto-spectra, freq x mic

if opt.legacy
    g2 = min(abs(Sxy) .^ 2 ./ (Sxx(:, i) .* Sxx(:, j) + eps), 0.99);
    Wt = g2 ./ (1 - g2);
    G = Wt .* Sxy ./ (abs(Sxy) + eps);
    snr_db = NaN;
else
    C = Sxy;
    useNoise = ~isempty(noise) && isequal(size(noise.Pxy), size(Sxy));
    if useNoise
        C = Sxy - noise.Pxy;                          % remove the stationary part
        snr_db = 10 * log10(mean(Sxx(:)) / max(mean(noise.Pxx(:)), eps));
    else
        snr_db = NaN;
    end
    g2 = abs(C) .^ 2 ./ (Sxx(:, i) .* Sxx(:, j) + eps) - 1.1 / nfr;   % bias-corrected
    g2 = min(max(g2, 0), 0.99);
    Wt = g2 ./ (1 - g2);
    G = Wt .* C ./ (abs(C) + eps);
end

step = 1 / (fs * opt.up);
Lmax = max(maxlag) + 2 / fs;
lags = (-ceil(Lmax / step):ceil(Lmax / step)) * step;
E = exp(1j * 2 * pi * f * lags);                      % freq x lag
R = zeros(numel(lags), P);
wsum = zeros(1, P);
for k = 1:size(opt.bands, 1)
    b = f >= opt.bands(k, 1) & f < opt.bands(k, 2);
    if ~any(b), continue; end
    Rk = real(E(b, :).' * G(b, :));                   % lag x pair
    Rk = Rk ./ max(max(abs(Rk), [], 1), eps);
    if opt.legacy
        ek = ones(1, P);
    else
        ek = sum(Wt(b, :), 1);                        % evidence of this band, per pair
    end
    R = R + Rk .* ek;
    wsum = wsum + ek;
end
R = R ./ max(wsum, eps);

info = struct('f', f, 'Sxx', Sxx, 'Sxy', Sxy, 'nfr', nfr, 'R', R, 'lags', lags, ...
    'evid', sum(Wt, 1), 'snr_db', snr_db);

Rm = R;
Rm(abs(lags(:)) > maxlag(:)' + 1 / fs) = -inf;        % physically impossible lags
[pk, idx] = max(Rm, [], 1);
tau = lags(idx)';
for p = 1:numel(idx)
    k = idx(p);
    if k > 1 && k < numel(lags) && all(isfinite(Rm(k - 1:k + 1, p)))
        a = Rm(k - 1, p); b = Rm(k, p); c = Rm(k + 1, p);
        den = a - 2 * b + c;
        if abs(den) > eps, tau(p) = tau(p) + 0.5 * (a - c) / den * step; end
    end
end
conf = max(pk(:), 0) .* mean(g2, 1)';
end
