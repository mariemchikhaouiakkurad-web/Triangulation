function out = tdoa_srp(info, W, B, grid)
% TDOA_SRP  Steered response: ALL pairs of ALL boards vote for every candidate source point
% before any decision (no per-pair argmax first).  Pairs stay inside one board, so the
% boards still need no clock sync.  This works at much lower SNR than picking 142 separate
% TDOA peaks and fitting them: a weak sound gives no clear peak on any single pair, but the
% sum over pairs still peaks at the right direction.
%
%   info  1 x nb cell, third output of tdoa_estimate per board ([] = board skipped)
%   W     1 x nb cell, world mic positions (mics x 3, m) from tdoa_pose
%   B     tdoa_boards(cfg) (pairs, weight)
%   grid  .az .el (deg, az 0 = +Y, +90 = +X)   .r  ranges (m, Inf = far field)
%         .fine_az .fine_el  refinement steps (deg) around the peak
%   out   .az .el .r          best point (r only meaningful below ~1 m with this aperture)
%         .peak               SRP value at the peak (consensus, about -1 .. 1)
%         .z                  peak contrast: (peak - median) / robust spread of the map
%         .map                az x el map at the best range
%         .wb                 weight each board got (board_weight x its SNR evidence)
c = 343;
nb = numel(B);
[AZ, EL] = ndgrid(grid.az, grid.el);
U = dirs(AZ(:), EL(:));
nr = numel(grid.r);

wb = zeros(1, nb);
for b = 1:nb
    if isempty(info{b}), continue; end
    m = mean(info{b}.evid) / max(numel(info{b}.f), 1);   % mean weight per bin = SNR-like
    wb(b) = B(b).weight * m / (1 + m);
end
if sum(wb) <= 0
    out = struct('az', NaN, 'el', NaN, 'r', NaN, 'peak', 0, 'z', 0, ...
        'map', zeros(numel(grid.az), numel(grid.el)), 'wb', wb);
    return
end

V = zeros(size(U, 1), nr);
for ir = 1:nr
    V(:, ir) = srp_at(U, grid.r(ir), info, W, B, wb, c);
end
[pk, k] = max(V(:));
[ku, ir] = ind2sub(size(V), k);
v = V(:, ir);
med = median(v);
spread = 1.4826 * median(abs(v - med)) + eps;

% refine around the peak at the best range
az0 = AZ(ku); el0 = EL(ku);
daz = abs(grid.az(min(2, end)) - grid.az(1)); del = abs(grid.el(min(2, end)) - grid.el(1));
[FA, FE] = ndgrid(az0 + (-daz:grid.fine_az:daz), ...
                  max(min(grid.el), el0 - del):grid.fine_el:min(max(grid.el), el0 + del));
Uf = dirs(FA(:), FE(:));
vf = srp_at(Uf, grid.r(ir), info, W, B, wb, c);
[pkf, kf] = max(vf);
if pkf >= pk
    az0 = FA(kf); el0 = FE(kf); pk = pkf;
end
az0 = mod(az0 + 180, 360) - 180;

out = struct('az', az0, 'el', el0, 'r', grid.r(ir), 'peak', pk, ...
    'z', (pk - med) / spread, 'map', reshape(v, numel(grid.az), numel(grid.el)), 'wb', wb);
end

function U = dirs(az, el)
U = [sind(az) .* cosd(el), cosd(az) .* cosd(el), sind(el)];
end

function v = srp_at(U, r, info, W, B, wb, c)
% weighted mean over boards of the mean over pairs of each pair's curve at its predicted lag
K = size(U, 1);
v = zeros(K, 1);
for b = 1:numel(B)
    if wb(b) <= 0, continue; end
    Wb = W{b}; pr = B(b).pairs; P = size(pr, 1);
    if isinf(r)
        tau = -(U * (Wb(pr(:, 1), :) - Wb(pr(:, 2), :))') / c;           % far field, K x P
    else
        Ps = r * U;                                                       % source points
        D = sqrt(max(sum(Ps .^ 2, 2) + sum(Wb .^ 2, 2)' - 2 * Ps * Wb', 0));   % K x mics
        tau = (D(:, pr(:, 1)) - D(:, pr(:, 2))) / c;
    end
    lags = info{b}.lags; L = numel(lags); step = lags(2) - lags(1);
    x = (tau - lags(1)) / step + 1;
    x = min(max(x, 1), L - 1e-9);
    k0 = floor(x); a = x - k0;
    idx = k0 + (0:P - 1) * L;
    Rb = info{b}.R;
    v = v + wb(b) * mean((1 - a) .* Rb(idx) + a .* Rb(idx + 1), 2);
end
v = v / sum(wb);
end
