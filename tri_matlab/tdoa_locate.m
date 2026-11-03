function [p, info] = tdoa_locate(B, tau, conf, cfg)
% TDOA_LOCATE  Stack the pairs of all calibrated boards and solve for the source position.
%   B        boards from tdoa_boards with .pose set by tdoa_calibrate
%   tau/conf 1 x 3 cells, per board TDOAs (s) and confidences
%   p        [x y z] (m);  info.std_xy  [std_x std_y std_z] (m), std_z = 0 in 2-D mode
%            info.cov  3 x 3 covariance (m^2), info.k95  scale of the 95 % ellipsoid (2.80 sigma)
%            info.board_us  per-board consistency residual (us), info.quality  weight factors
%
% Board quality (cfg.tdoa.adaptive): the TDOAs of one board must fit ONE direction.  The
% plane-wave residual is computed in world coordinates, so it does not depend on how the
% board is oriented.  A close source bends the wavefront and adds a residual that grows with
% the board size; that part is predicted at the solved position and removed, so a large
% board is not penalised for its size.  A board whose remaining residual is large (echoes,
% obstruction) is down-weighted by (best residual / its residual)^2 relative to the best
% board of this block, then the position is solved again.
n = numel(B);
W = cell(1, n); D = cell(1, n); res_us = nan(1, n);
for b = 1:n
    W{b} = tdoa_pose(B(b).Mloc, B(b).pose);
    D{b} = W{b}(B(b).pairs(:, 1), :) - W{b}(B(b).pairs(:, 2), :);
    if size(W{b}, 1) >= 3, res_us(b) = planefit_us(D{b}, tau{b}(:)); end
end
adaptive = isfield(cfg.tdoa, 'adaptive') && cfg.tdoa.adaptive;
q = quality(res_us, adaptive);
[p, info, C] = solve(B, W, tau, conf, cfg, q);
if adaptive && any(isfinite(res_us))
    % remove the curvature part predicted at the solved position, then re-weight
    c = 343; cur_us = zeros(1, n);
    for b = find(isfinite(res_us))
        d = vecnorm(p - W{b}, 2, 2);
        cur_us(b) = planefit_us(D{b}, (d(B(b).pairs(:, 1)) - d(B(b).pairs(:, 2))) / c);
    end
    k = isfinite(res_us);                             % boards without a residual (2 mics) stay NaN:
    res_us(k) = sqrt(max(res_us(k) .^ 2 - cur_us(k) .^ 2, 0));   % max(NaN, 0) would give 0 = "perfect"
    q2 = quality(res_us, adaptive);
    if max(abs(q2 - q)) > 0.1
        q = q2;
        [p, info, C] = solve(B, W, tau, conf, cfg, q);
    end
end
info.cov = C;
info.board_us = res_us;
info.quality = q;
info.std_xy = sqrt(max(diag(C), 0))';
info.k95 = 2.7955;                                    % sqrt(chi2inv(0.95, 3)): 95 % of the true positions
end

% =====================================================================================
function r = planefit_us(D, tau)
% RMS misfit (us) of the best single plane-wave direction to one board's TDOAs
s = -pinv(D) * tau;
r = rms(tau + D * s) * 1e6;
end

function q = quality(res_us, adaptive)
q = ones(size(res_us));
if adaptive && any(isfinite(res_us))
    ref = max(min(res_us), 5);                        % best board of this block (>= 5 us)
    k = isfinite(res_us);
    q(k) = (ref ./ max(res_us(k), ref)) .^ 2;
end
end

function [p, info, C] = solve(B, W, tau, conf, cfg, q)
Mi = []; Mj = []; tt = []; ww = [];
infl = 0;
for b = 1:numel(B)
    c = conf{b}(:);
    Mi = [Mi; W{b}(B(b).pairs(:, 1), :)]; %#ok<AGROW>
    Mj = [Mj; W{b}(B(b).pairs(:, 2), :)]; %#ok<AGROW>
    tt = [tt; tau{b}(:)]; %#ok<AGROW>
    ww = [ww; B(b).weight * q(b) * c / max(sum(c), eps)]; %#ok<AGROW>   % board total = weight x quality
    % n mics give n(n-1)/2 pairs but only n-1 independent TDOAs: counting all pairs as
    % independent overstates the information about n/2 times
    infl = infl + B(b).weight * q(b) * size(B(b).Mloc, 1) / 2;
end
infl = infl / sum([B.weight] .* q);
opt = cfg.tdoa;
opt.grid_x = cfg.tri.grid_x; opt.grid_y = cfg.tri.grid_y;
opt.exclude = vertcat(cfg.tri.pos{:});
opt.origin = mean(opt.exclude, 1);                    % centre of the boards (far-field search)
if numel(opt.origin) == 2, opt.origin(3) = 0; end
[p, info] = tdoa_solve(tt, ww, Mi, Mj, opt);
C = info.cov * infl;
if numel(p) == 2                                      % 2-D mode: source at the fixed height src_z
    p = [p, opt.src_z];
    C = blkdiag(C, 0);
end
end
