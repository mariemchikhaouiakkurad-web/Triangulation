function dir = freq_direction(spec, det, B, F)
% FREQ_DIRECTION  Direction of the detected harmonic comb, from its harmonics >= dir_fmin only.
%
%   spec  1 x nb cell, bins x mics x frames spectra (from freq_detect)
%   det   detection (freq_detect): harm_f, harm_bins (bin per frame), harm_snr
%   B     boards (tdoa_boards) with .pose from the calibration (world mic positions via tdoa_pose)
%   F     cfg.freq: dir_fmin (Hz), el_range [min max] (deg), az_step (deg)
%   dir   .az .el (deg; az 0 = +Y toward the UMA-8, +90 = +X right; el above the boards),
%         NaN when no harmonic can be used (.reason says why)
%         .daz .del  half-width of the peak at -1 dB (deg)   .ambiguity  main / 2nd peak (dB)
%         .ambiguous 2nd peak within 1.5 dB                 .map  az x el power map (0..1)
%         .used      harmonics used (logical, per harmonic) .elec  harmonics found electrical
%
% Per board b and harmonic k: cross-spectral matrix R = mean over frames of x x^H (x taken at the
% harmonic's bin in each frame, so a drifting comb is followed), normalised element by element
% (PHAT).  Steered power P_bk(u) = v^H R v / M^2, v_m = exp(+j 2 pi f (w_m . u) / c), w_m = the
% calibrated world mic positions of that board.  P = sum over boards and harmonics, weighted by
% the harmonic's SNR and the board weight.  Each board only uses its own mics: no sync needed.
% Harmonics below dir_fmin are not used: their wavelength (> 0.4 m) is much larger than the boards,
% so they barely tell directions apart (they still count for the detection).
%
% Electrical test: an electrical tone (USB, power supply, laptop) reaches every mic with no delay,
% so every board sees zero phase between its mics at once.  A real sound can make one board look
% like that (source straight in front of it) but not all boards together, since they face
% different ways.  A harmonic is electrical when the zero-delay fit of the multi-mic boards beats
% the best real direction by > 0.05 (and is > 0.9); it is then left out of the direction.
c = 343;
azg = -180:F.az_step:180 - F.az_step;
elg = F.el_range(1):F.az_step:F.el_range(2);
[AZ, EL] = meshgrid(azg, elg);
nh = numel(det.harm_f);
W = cell(1, numel(B)); for b = 1:numel(B), W{b} = tdoa_pose(B(b).Mloc, B(b).pose); end
R = cell(numel(B), nh);
for b = 1:numel(B)
    for k = 1:nh
        Xk = zeros(size(W{b}, 1), size(det.harm_bins, 2));
        for t = 1:size(det.harm_bins, 2), Xk(:, t) = spec{b}(det.harm_bins(k, t), :, t).'; end
        Rk = (Xk * Xk') / size(Xk, 2);
        R{b, k} = Rk ./ max(abs(Rk), eps);                     % PHAT: phases only
    end
end
% electrical test, per harmonic, with the multi-mic boards
multi = find(arrayfun(@(b) size(W{b}, 1) >= 3, 1:numel(B)));
dir.elec = false(1, nh);
for k = 1:nh
    zp = 0; best = zeros(numel(AZ), 1); ws = 0;
    for b = multi
        M = size(W{b}, 1);
        zp = zp + B(b).weight * real(sum(R{b, k}(:))) / M ^ 2;      % steering with zero delays
        best = best + B(b).weight * steer1(W{b}, R{b, k}, det.harm_f(k), AZ(:), EL(:));
        ws = ws + B(b).weight;
    end
    zp = zp / ws; best = max(best) / ws;
    dir.elec(k) = zp > 0.9 && zp - best > 0.05;
end
dir.used = det.harm_f >= F.dir_fmin & ~dir.elec;
dir.map = []; dir.az_grid = azg; dir.el_grid = elg;
dir.az = NaN; dir.el = NaN; dir.daz = NaN; dir.del = NaN; dir.ambiguity = NaN; dir.ambiguous = true;
if ~any(dir.used)
    if any(dir.elec), dir.reason = 'harmonics electrical';
    else, dir.reason = sprintf('no harmonic >= %d Hz', F.dir_fmin); end
    return
end
dir.reason = '';
P = steer(AZ(:), EL(:));
P = P / max(P);
dir.map = reshape(P, size(AZ));
[pmax, i] = max(P);
% refine around the coarse peak
[A2, E2] = meshgrid(AZ(i) + (-F.az_step:0.5:F.az_step), ...
    min(max(EL(i) + (-F.az_step:0.5:F.az_step), -90), 90));
P2 = steer(A2(:), E2(:));
[~, j] = max(P2);
dir.az = mod(A2(j) + 180, 360) - 180; dir.el = E2(j);
% width of the main peak (-1 dB) on the coarse map
dist = angdist(AZ(:), EL(:), AZ(i), EL(i));
in1 = P >= pmax * 10 ^ (-1 / 10) & dist <= 90;
daz = abs(mod(AZ(in1) - AZ(i) + 180, 360) - 180);
dir.daz = max(daz) + F.az_step / 2;
dir.del = max(abs(EL(in1) - EL(i))) + F.az_step / 2;
% strongest OTHER peak: outside the main lobe (where the power has dropped to half)
lobe = max(dist(P >= pmax / 2 & dist <= 120));
p2 = max([P(dist > max(lobe, 15)); eps]);
dir.ambiguity = 10 * log10(pmax / p2);
dir.ambiguous = dir.ambiguity < 1.5;

    function P = steer(az, el)
        P = zeros(numel(az), 1);
        for bb = 1:numel(B)
            for kk = find(dir.used)
                P = P + B(bb).weight * max(det.harm_snr(kk), 0) * steer1(W{bb}, R{bb, kk}, det.harm_f(kk), az, el);
            end
        end
    end
end

function p = steer1(Wb, R, f, az, el)
% normalised steered power (0..1) of one board at frequency f toward the directions az, el
c = 343;
u = [cosd(el) .* sind(az), cosd(el) .* cosd(az), sind(el)];
V = exp(1j * 2 * pi * f * (u * Wb') / c);
p = real(sum(conj(V) .* (V * R.'), 2)) / size(Wb, 1) ^ 2;
end

function d = angdist(az1, el1, az2, el2)
% great-circle angle between directions (deg)
d = acosd(min(1, max(-1, sind(el1) .* sind(el2) + cosd(el1) .* cosd(el2) .* cosd(az1 - az2))));
end
