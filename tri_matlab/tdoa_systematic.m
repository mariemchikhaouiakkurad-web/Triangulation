function sys = tdoa_systematic(B, taucal, cal, cfg)
% TDOA_SYSTEMATIC  Systematic error of the method, measured at the calibration points.
%
%   taucal  {point, board} median TDOAs at every calibration point
%   cal     n x 3 true positions of the calibration points (m)
%   sys     .range  RMS error ALONG the line of sight (distance error, m)
%           .perp   RMS error ACROSS the line of sight, per axis (direction error, m)
%           .bias   mean distance error (m), > 0 = measured too far, < 0 = too close
%           .centre the point the line of sight starts from (centre of the boards)
%
% The statistical uncertainty of tdoa_solve only covers the noise of the TDOAs; errors of
% the geometry (board positions, heights, calibration) add a bias it cannot see.  This
% bias is measured where the true position is known and added to the 95 % region.
sys.centre = mean(vertcat(cfg.tri.pos{:}), 1);
if numel(sys.centre) == 2, sys.centre(3) = 0; end
n = size(cal, 1);
er = zeros(n, 1); ep = zeros(n, 1);
fprintf('\nerror at the calibration points (known positions):\n');
for p = 1:n
    conf = cellfun(@(t) ones(size(t)), taucal(p, :), 'UniformOutput', false);
    q = tdoa_locate(B, taucal(p, :), conf, cfg);
    u = cal(p, :) - sys.centre; u = u / norm(u);
    e = q - cal(p, :);
    er(p) = e * u';                                   % along the line of sight
    ep(p) = norm(e - er(p) * u);                      % across it
    fprintf('  cal %d  true (%+.2f %+.2f %+.2f)  measured (%+.2f %+.2f %+.2f)  distance %+4.0f cm, direction %3.0f cm\n', ...
        p, cal(p, :), q, 100 * er(p), 100 * ep(p));
end
% robust (median) so one badly placed or badly heard point does not blow up the 95 % region
sys.range = sqrt(median(er .^ 2));
sys.perp = sqrt(median(ep .^ 2) / 2);
sys.bias = median(er);
bad = find(hypot(er, ep) > 3 * max(hypot(sys.range, sqrt(2) * sys.perp), 0.05));
fprintf('systematic error: distance %.0f cm (median, bias %+.0f cm), direction %.0f cm per axis\n', ...
    100 * sys.range, 100 * sys.bias, 100 * sys.perp);
if ~isempty(bad)
    fprintf('  outlier calibration point(s): %s - phone not at the point, or poorly covered by the boards\n', ...
        strjoin(compose('cal %d', bad), ', '));
end
end
