function [az, el, score] = tri_srp(X, Mpos, fs, cfg)
% TRI_SRP  SRP-PHAT direction of one block X (samples x mics, positions Mpos in m),
% coarse grid then fine around the peak (as in tri_tdoa_combine), plus a 0-1 score.
% az 0 = +y, +90 = +x; el = angle above the x-y plane.
nf = cfg.nfft; hop = nf / 2;
w = 0.5 - 0.5 * cos(2 * pi * (0:nf - 1)' / nf);
f = (0:nf / 2)' * fs / nf;
sel = find(f >= cfg.band(1) & f <= cfg.band(2));
fr = f(sel);
nfr = floor((size(X, 1) - nf) / hop) + 1;
Y = zeros(numel(sel), size(X, 2), nfr);
for k = 1:nfr
    F = fft(X((k - 1) * hop + (1:nf), :) .* w);
    F = F(sel, :);
    Y(:, :, k) = F ./ (abs(F) + 1e-12);
end
[AZ, EL] = meshgrid(cfg.grid_az, cfg.grid_el);
AZ = AZ(:); EL = EL(:);
if isfield(cfg, 'upper_only') && cfg.upper_only
    AZ = AZ(EL >= 0); EL = EL(EL >= 0);
end
P = srp_eval(Y, Mpos, fr, AZ, EL, cfg);
[~, i] = max(P);
[AZ, EL] = meshgrid(AZ(i) + (-4:0.5:4), max(-90, min(90, EL(i) + (-6:0.5:6))));
P = srp_eval(Y, Mpos, fr, AZ(:), EL(:), cfg);
[pmax, i] = max(P);
az = mod(AZ(i) + 180, 360) - 180; el = EL(i);
score = pmax / (numel(fr) * nfr * size(Mpos, 1) ^ 2);
end

function P = srp_eval(Y, Mpos, fr, az, el, cfg)
c0 = 343;
u = [sind(az) .* cosd(el), cosd(az) .* cosd(el), sind(el)];
if isfield(cfg, 'range') && isfinite(cfg.range)
    S = mean(Mpos, 1) + cfg.range * u;
    tau = zeros(size(u, 1), size(Mpos, 1));
    for m = 1:size(Mpos, 1)
        tau(:, m) = sqrt(sum((S - Mpos(m, :)) .^ 2, 2)) / c0;
    end
else
    tau = -(u * Mpos') / c0;
end
P = zeros(size(u, 1), 1);
for fi = 1:numel(fr)
    A = exp(1j * 2 * pi * fr(fi) * tau);
    Yf = reshape(Y(fi, :, :), size(Y, 2), []);
    P = P + sum(abs(A * Yf) .^ 2, 2);
end
end
