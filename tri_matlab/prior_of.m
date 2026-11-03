function p = prior_of(S, b)
% PRIOR_OF  Mounting prior of board b from cfg.tdoa (S), or [] when none is configured.
if isfield(S, 'prior') && numel(S.prior) >= b
    p = S.prior(b);
else
    p = [];
end
end
