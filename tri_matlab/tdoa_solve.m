function [p, info] = tdoa_solve(tau, w, Mi, Mj, opt)
% TDOA_SOLVE  Source position from time differences (hyperbolic multilateration).
%
%   tau     P x 1 measured tau = arrival at mic i - arrival at mic j (s)
%   w       P x 1 weights (confidence)
%   Mi, Mj  P x 3 world positions (m) of the two mics of every pair
%   opt     .grid_x, .grid_y  search area (m)
%           .solve_z  true -> 3-D (x, y, z) inside .grid_z;  false -> 2-D at height .src_z
%           .max_range  (3-D) search directions x distances 0.25..max_range m around .origin,
%                       elevations .el_grid (deg), instead of the x-y-z box
%           .coarse  start-grid step (m) in 3-D
%           .sigma  expected TDOA error (s)       .huber  robust threshold (in sigmas)
%           .iters  Levenberg-Marquardt iterations
%           .exclude  K x 2 points to avoid (boards, 10 cm radius in x-y)
%   p       [x y] or [x y z] position (m)
%   info    .std_xy  1-sigma uncertainty per coordinate (m)   .cov  covariance (m^2)
%           .rms_us  weighted residual RMS (us)   .resid  P x 1 residuals (s)
%
% Model:  tau_ij(p) = (|p - m_i| - |p - m_j|) / c   (near field, exact)
%   1. coarse grid: robust (Huber) cost at every grid point -> start point
%   2. Levenberg-Marquardt with iteratively re-weighted Huber weights -> refined position
%   3. covariance (J' W J)^-1 scaled by the residual variance -> uncertainty
c = 343;
w = w(:) / max(sum(w), eps);
solve_z = isfield(opt, 'solve_z') && opt.solve_z;
if solve_z && isfield(opt, 'max_range') && ~isempty(opt.max_range)
    % far field: start grid in direction x distance around opt.origin, out to max_range
    [az, el, rr] = ndgrid(-180:4:176, opt.el_grid, ...
        logspace(log10(0.25), log10(opt.max_range), 24));
    Q = opt.origin + rr(:) .* [cosd(el(:)) .* sind(az(:)), cosd(el(:)) .* cosd(az(:)), sind(el(:))];
    lo = opt.origin - opt.max_range; hi = opt.origin + opt.max_range;
    lo(3) = opt.origin(3) - 0.10;                     % sources are not under the table
    hi(3) = opt.origin(3) + opt.max_range * sind(opt.el_grid(end)) + 0.05;
    Q = Q(Q(:, 3) >= lo(3), :);                       % start points inside the height limits too
elseif solve_z
    s = opt.coarse;
    [gx, gy, gz] = ndgrid(opt.grid_x(1):s:opt.grid_x(end), opt.grid_y(1):s:opt.grid_y(end), ...
                          opt.grid_z(1):s:opt.grid_z(end));
    Q = [gx(:), gy(:), gz(:)];
    lo = [opt.grid_x(1), opt.grid_y(1), opt.grid_z(1)];
    hi = [opt.grid_x(end), opt.grid_y(end), opt.grid_z(end)];
else
    [gx, gy] = meshgrid(opt.grid_x, opt.grid_y);
    Q = [gx(:), gy(:)];
    lo = [opt.grid_x(1), opt.grid_y(1)]; hi = [opt.grid_x(end), opt.grid_y(end)];
end
keep = true(size(Q, 1), 1);
for k = 1:size(opt.exclude, 1)
    keep = keep & hypot(Q(:, 1) - opt.exclude(k, 1), Q(:, 2) - opt.exclude(k, 2)) >= 0.10;
end
Q = Q(keep, :);

% 1. coarse grid
cost = zeros(size(Q, 1), 1);
for s = 1:4096:size(Q, 1)
    q = Q(s:min(end, s + 4095), :);
    r = (model(q, Mi, Mj, opt.src_z, c) - tau') / opt.sigma;   % points x pairs
    cost(s:s + size(q, 1) - 1) = huber(r, opt.huber) * w;
end
[~, i0] = min(cost);
p = Q(i0, :);

% 2. Levenberg-Marquardt + IRLS: a step is kept only if the robust cost goes down and
%    the point stays inside the search area; otherwise the damping is increased.
f0 = robust_cost(p, tau, w, Mi, Mj, opt, c);
lambda = 1e-2;
for it = 1:opt.iters
    [pred, J] = model(p, Mi, Mj, opt.src_z, c);
    r = tau - pred';
    W = w .* min(1, opt.huber * opt.sigma ./ max(abs(r), eps));   % Huber weights
    A = J' * (W .* J);
    g = J' * (W .* r);
    improved = false;
    while lambda < 1e6
        dp = ((A + lambda * diag(diag(A) + eps)) \ g)';
        q = min(max(p + dp, lo), hi);
        f1 = robust_cost(q, tau, w, Mi, Mj, opt, c);
        if f1 < f0
            p = q; f0 = f1; lambda = max(lambda / 3, 1e-6); improved = true; break;
        end
        lambda = lambda * 4;
    end
    if ~improved || norm(dp) < 1e-4, break; end
end

if isfield(opt, 'max_range') && ~isempty(opt.max_range) && numel(p) == 3
    d = p - opt.origin;                               % keep inside the search sphere, not its box
    if norm(d) > opt.max_range, p = opt.origin + d * opt.max_range / norm(d); end
end

% 3. uncertainty
[pred, J] = model(p, Mi, Mj, opt.src_z, c);
r = tau - pred';
W = w .* min(1, opt.huber * opt.sigma ./ max(abs(r), eps));
s2 = sum(W .* r .^ 2) / max(sum(W), eps);
C = s2 * pinv(J' * (W .* J)) * sum(W .^ 2) / max(sum(W) ^ 2, eps);
info.cov = (C + C') / 2;
info.std_xy = sqrt(max(diag(C), 0))';
info.rms_us = sqrt(s2) * 1e6;
info.resid = r;
end

function [t, J] = model(q, Mi, Mj, z, c)
% predicted TDOAs for points q (K x 2, at height z, or K x 3) and, for one point,
% the Jacobian d tau / d(q)
nq = size(q, 2);
if nq == 2, q3 = [q, z * ones(size(q, 1), 1)]; else, q3 = q; end
di = sqrt((q3(:, 1) - Mi(:, 1)') .^ 2 + (q3(:, 2) - Mi(:, 2)') .^ 2 + (q3(:, 3) - Mi(:, 3)') .^ 2);
dj = sqrt((q3(:, 1) - Mj(:, 1)') .^ 2 + (q3(:, 2) - Mj(:, 2)') .^ 2 + (q3(:, 3) - Mj(:, 3)') .^ 2);
t = (di - dj) / c;
if nargout > 1
    J = ((q3(1:nq) - Mi(:, 1:nq)) ./ di' - (q3(1:nq) - Mj(:, 1:nq)) ./ dj') / c;
end
end

function f = robust_cost(p, tau, w, Mi, Mj, opt, c)
f = huber((model(p, Mi, Mj, opt.src_z, c) - tau') / opt.sigma, opt.huber) * w;
end

function rho = huber(r, k)
a = abs(r);
rho = (a <= k) .* (0.5 * a .^ 2) + (a > k) .* (k * a - 0.5 * k ^ 2);
end
