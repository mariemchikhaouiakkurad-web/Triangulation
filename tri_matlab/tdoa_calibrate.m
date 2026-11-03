function [pose, err_us] = tdoa_calibrate(Mloc, centre, pairs, taucal, calpts, src_z, flat, prior)
% TDOA_CALIBRATE  3-D orientation (+ scale) of one board from TDOAs measured with the
% phone at known points.  No ruler measurement of the board is needed.
%
%   Mloc    mics x 3 board-frame positions (m), planar board (z = 0)
%   centre  [x y] or [x y z] rough board centre in the world (m)
%   pairs   P x 2 mic pairs;  taucal  P x n measured TDOAs (s), one column per point
%   calpts  n x 2 (or n x 3) known phone positions (m);  src_z  phone height if n x 2
%   flat    true = the board is known to lie flat (only rotation about z, mirror, scale)
%   prior   (optional) what is known about the mounting:
%           .normal  1x3 direction the mic face points to (front or back: the sign is free,
%                    a planar board cannot tell them apart);  [] = unknown
%           .tol     the fitted normal must stay within tol degrees of it
%           .scale   fixed geometry scale (1 = layout exact);  NaN = fitted
%           For a 2-mic board the MIC1-MIC16 axis is kept in the board plane (normal within tol).
%   pose    .R  3x3 rotation board -> world   .scale   .centre   (used by tdoa_pose)
%   err_us  RMS misfit of the TDOAs at the calibration points (us)
%
% Far field: tau_ij = -(w_i - w_j) . u / c, w = world mic positions, u = unit vector toward
% the phone.  For a planar board only the two in-plane axes r1, r2 of R matter, so per point
% the in-plane slowness s = -D \ tau (D = in-plane baselines) must equal scale * [r1 r2]' u / c.
%   1. coarse search over yaw / pitch / roll (and mirror: r1 -> -r1), scale in closed form
%   2. fminsearch refinement
%   3. near-field refinement (rotation, scale, board centre within ~3 cm) with the exact model
%      (not for the 2-mic board, and the origin board's centre stays at 0)
%   4. misfit with the exact near-field model at the calibration points
% With a prior, the search starts only near the known facing and is penalised beyond tol.
% Use >= 3 points: with 2 points a tilted board is only just determined.
c = 343;
if nargin < 7, flat = false; end
if nargin < 8 || isempty(prior), prior = struct('normal', [], 'tol', 90, 'scale', NaN); end
nrm = prior.normal; if ~isempty(nrm), nrm = nrm(:) / norm(nrm); end
if numel(centre) == 2, centre = [centre(:)', 0]; end
if size(calpts, 2) == 2, calpts = [calpts, src_z * ones(size(calpts, 1), 1)]; end
n = size(calpts, 1);
U = (calpts - centre)'; U = U ./ sqrt(sum(U .^ 2, 1));            % 3 x n directions

D = Mloc(pairs(:, 1), 1:2) - Mloc(pairs(:, 2), 1:2);              % in-plane baselines
twoMics = size(Mloc, 1) == 2;
if twoMics
    % one pair: only its axis matters.  tau_k = -(scale*d*a) . u_k / c, a = unit axis in the world
    [az, el] = meshgrid(-180:2:178, -90:2:90);
    A = [cosd(el(:)) .* sind(az(:)), cosd(el(:)) .* cosd(az(:)), sind(el(:))];   % candidate axes
    if ~isempty(nrm)                                              % the axis lies in the board plane
        A = A(abs(A * nrm) <= sind(prior.tol), :);
    end
    d = norm(D);
    pred = -(d / c) * (A * U);                                    % candidates x n
    e = sum((pred - taucal(1, :)) .^ 2, 2);
    [~, i] = min(e);
    ax = A(i, :)';
    % board frame x axis = MIC1 -> MIC16 direction; build any rotation with r1 along the axis
    lx = D(1, :)' / d; lx = [lx; 0];
    r1 = ax; r3 = null(r1'); r3 = r3(:, 1); r2 = cross(r3, r1);
    Rw = [r1 r2 r3]; Rl = [lx, [-lx(2); lx(1); 0], [0; 0; 1]];
    pose = struct('R', Rw * Rl', 'scale', 1, 'centre', centre);
else
    S = -(D \ taucal);                                            % 2 x n measured slowness
    fix = prior.scale;
    if flat                                                       % yaw only: pitch = roll = 0
        cost = @(q, m) slowness_misfit([q(1) 0 0], m, S, U, c, fix, nrm, prior.tol);
        tilts = 0;
    else
        cost = @(q, m) slowness_misfit(q, m, S, U, c, fix, nrm, prior.tol);
        tilts = -90:15:90;
    end
    best = inf;
    for m = [1 -1]
        [Y, Pp, Rr] = ndgrid(0:10:350, tilts, tilts);
        for k = 1:numel(Y)
            if ~isempty(nrm)                                      % start only near the known facing
                R0 = rotm([Y(k) Pp(k) Rr(k)]);
                if acosd(min(1, abs(R0(:, 3)' * nrm))) > prior.tol + 15, continue; end
            end
            e = cost([Y(k) Pp(k) Rr(k)], m);
            if e < best, best = e; q0 = [Y(k) Pp(k) Rr(k)]; m0 = m; end
        end
    end
    if flat, q0 = q0(1); end
    q = fminsearch(@(q) cost(q, m0), q0, optimset('TolX', 1e-3, 'TolFun', 1e-14, 'Display', 'off'));
    if flat, q = [q 0 0]; end
    [~, k] = cost(q, m0);
    R = rotm(q); R(:, 1) = m0 * R(:, 1);
    pose = struct('R', R, 'scale', k, 'centre', centre);
    % 4. near-field refinement: the far-field fit above is biased at ~1 m (wavefront sagitta
    %    across a 13 cm board ~2 mm = ~6 us, same order as cfg.tdoa.sigma).  Re-fit rotation,
    %    scale and the board centre (within ~3 cm of the ruler value) with the exact model.
    pose = nf_refine(Mloc, pose, q, m0, k, fix, flat, centre, pairs, taucal, calpts, nrm, prior.tol, c);
end

% misfit with the exact near-field model
W = tdoa_pose(Mloc, pose);
e = 0;
for k = 1:n
    dd = sqrt(sum((W - calpts(k, :)) .^ 2, 2));
    e = e + sum(((dd(pairs(:, 1)) - dd(pairs(:, 2))) / c - taucal(:, k)) .^ 2);
end
err_us = sqrt(e / numel(taucal)) * 1e6;
end

function [e, k] = slowness_misfit(q, m, S, U, c, fix, nrm, tol)
R = rotm(q);
V = [m * R(:, 1), R(:, 2)]' * U / c;                             % 2 x n, before scale
if isfinite(fix)
    k = fix;                                                      % known geometry scale
else
    k = sum(S(:) .* V(:)) / max(sum(V(:) .^ 2), eps);           % best scale (closed form)
    k = min(max(k, 0.5), 2);
end
e = sum((S(:) - k * V(:)) .^ 2);
if ~isempty(nrm)                                                  % stay within tol of the known facing
    ang = acosd(min(1, abs(R(:, 3)' * nrm)));
    e = e + sum(S(:) .^ 2) * max(0, (ang - tol) / 5) ^ 2;
end
end

function pose = nf_refine(Mloc, pose0, q, m0, k, fix, flat, centre, pairs, taucal, calpts, nrm, tol, c)
% near-field least squares from the far-field solution; kept only if it lowers the misfit.
% The origin board (centre = [0 0 0], the UMA-8) defines the frame, so its centre stays fixed.
moveC = norm(centre) > 1e-9;
freeK = ~isfinite(fix);
if flat, x0 = q(1); else, x0 = q(:)'; end
nq = numel(x0);
if freeK, x0 = [x0, 0]; end                       % 100 x log(scale / far-field scale)
if moveC, x0 = [x0, 0 0 0]; end                   % centre shift (cm)
sig_us = 10;                                      % TDOA error scale for the centre prior
sig_cm = 3;                                       % ruler accuracy of cfg.tri.pos
opts = optimset('TolX', 1e-3, 'TolFun', 1e-6, 'MaxFunEvals', 4000, 'MaxIter', 4000, 'Display', 'off');
e0 = cost(x0);
x = fminsearch(@cost, x0, opts);
if cost(x) < e0
    pose = unpack(x);
else
    pose = pose0;
end

    function p = unpack(x)
        qq = x(1:nq); if flat, qq = [qq 0 0]; end
        kk = k; ii = nq;
        if freeK, ii = ii + 1; kk = min(max(k * exp(x(ii) / 100), 0.5), 2); end
        cc = centre;
        if moveC, cc = centre + x(ii + 1:ii + 3) / 100; end
        RR = rotm(qq); RR(:, 1) = m0 * RR(:, 1);
        p = struct('R', RR, 'scale', kk, 'centre', cc);
    end

    function e = cost(x)
        p = unpack(x);
        Wm = tdoa_pose(Mloc, p);
        e = 0;
        for kp = 1:size(calpts, 1)
            dd = sqrt(sum((Wm - calpts(kp, :)) .^ 2, 2));
            e = e + sum(((dd(pairs(:, 1)) - dd(pairs(:, 2))) / c - taucal(:, kp)) .^ 2);
        end
        e = e * 1e12;                                 % us^2
        if moveC
            e = e + sig_us ^ 2 * sum(((p.centre - centre) * 100 / sig_cm) .^ 2);
        end
        if ~isempty(nrm)
            ang = acosd(min(1, abs(p.R(:, 3)' * nrm)));
            e = e + numel(taucal) * sig_us ^ 2 * max(0, (ang - tol) / 5) ^ 2;
        end
    end
end

function R = rotm(q)
% yaw about z, then pitch about x, then roll about y (degrees)
y = q(1); p = q(2); r = q(3);
Rz = [cosd(y) -sind(y) 0; sind(y) cosd(y) 0; 0 0 1];
Rx = [1 0 0; 0 cosd(p) -sind(p); 0 sind(p) cosd(p)];
Ry = [cosd(r) 0 sind(r); 0 1 0; -sind(r) 0 cosd(r)];
R = Rz * Rx * Ry;
end
