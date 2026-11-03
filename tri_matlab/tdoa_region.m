function [C, axes95, V] = tdoa_region(p, Cstat, sys)
% TDOA_REGION  95 % region of the true position: statistical covariance + systematic error.
%
%   p       [x y z] estimate (m);  Cstat  3 x 3 covariance from tdoa_locate (m^2)
%   sys     from tdoa_systematic (range / perp / centre), or [] for statistics only
%   C       total covariance (m^2)
%   axes95  half-axes of the 95 % ellipsoid (m), smallest first;  V  their directions (columns)
%
% The systematic part is an ellipsoid aligned with the line of sight from the boards:
% sys.range along it (distance error), sys.perp across it (direction error).
C = Cstat;
if ~isempty(sys)
    u = p - sys.centre; u = u(:) / max(norm(u), eps);
    C = C + sys.perp ^ 2 * (eye(3) - u * u') + sys.range ^ 2 * (u * u');
end
C = (C + C') / 2;
[V, D] = eig(C);
[d, i] = sort(max(diag(D), 0));
V = V(:, i);
axes95 = 2.7955 * sqrt(d)';                           % sqrt(chi2inv(0.95, 3))
end
