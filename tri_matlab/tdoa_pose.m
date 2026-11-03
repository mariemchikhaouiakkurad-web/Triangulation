function W = tdoa_pose(Mloc, pose)
% TDOA_POSE  World positions of a board's mics: scale, rotate (board -> world), move to centre.
W = pose.scale * (pose.R * Mloc')' + pose.centre;
end
