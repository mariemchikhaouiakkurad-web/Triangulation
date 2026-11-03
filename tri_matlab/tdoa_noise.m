function noise = tdoa_noise(noise, info, alpha)
% TDOA_NOISE  Running noise model of ONE board (auto- and cross-spectra), as in freq mode's
% noise floor.  Updated only on blocks where nothing was detected, so a target is never learnt.
%
%   noise = tdoa_noise([], info)          start from the first block
%   noise = tdoa_noise(noise, info, a)    update; a = forgetting factor per block
%                                         (the first blocks use 1/n, i.e. a plain mean)
%   info  = third output of tdoa_estimate (same nfft / band)
%
%   noise.Pxx  freq x mic     noise.Pxy  freq x pair     noise.n  blocks used
if nargin < 3 || isempty(alpha), alpha = 0.03; end
if isempty(noise) || ~isequal(size(noise.Pxx), size(info.Sxx)) || ~isequal(size(noise.Pxy), size(info.Sxy))
    noise = struct('Pxx', info.Sxx, 'Pxy', info.Sxy, 'n', 1);
    return
end
a = max(alpha, 1 / (noise.n + 1));
noise.Pxx = (1 - a) * noise.Pxx + a * info.Sxx;
noise.Pxy = (1 - a) * noise.Pxy + a * info.Sxy;
noise.n = noise.n + 1;
end
