function tdoa_print_pose(B, err_us)
% TDOA_PRINT_POSE  One line per board: calibrated orientation, scale and misfit.
R = B.pose.R;
if size(B.Mloc, 1) == 2
    a = R * (B.Mloc(2, :) - B.Mloc(1, :))'; a = a / norm(a);
    fprintf('%-7s MIC1->MIC16 axis: az %+6.1f deg, el %+5.1f deg', B.name, atan2d(a(1), a(2)), asind(a(3)));
else
    nrm = R(:, 3);
    fprintf('%-7s tilt %4.0f deg (0 = flat)  x-axis az %+6.1f deg  scale %.2f', B.name, ...
        acosd(abs(nrm(3))), atan2d(R(1, 1), R(2, 1)), B.pose.scale);
end
fprintf('   misfit %5.1f us', err_us);
if err_us > 20, fprintf('  <-- poor: phone not at the point / too quiet'); end
fprintf('\n');
end
