% TRI_CONFIG  Your layout (2026-11-02): UMA-8 at the top, UMA16A (intact) on the LEFT, UMA16B (cut)
% on the RIGHT, laptop at the bottom.  The modes 'live3' / 'tdoa_*' use cfg.tri.pos below;
% the old 'live' / 'offline' modes use the MEASUREMENTS block.  Then:  >> tri_config
%
% Frame (metres):  O = centre of the UMA-8 (on the table: z = 0 is the table top)
%                  X = to the right as seen from the laptop
%                  Y = away from the laptop, toward the UMA-8
%                  Z = up (height above the table)
clear cfg
cfg.fs      = 48000;
cfg.names   = {'UMA16A', 'UMA8', 'UMA16B'};              % first = time reference
cfg.files   = {'uma16a.wav', 'uma8.wav', 'uma16b.wav'};  % recorded at the same time
cfg.odascfg = {'uma16.cfg', 'uma8.cfg', 'uma16.cfg'};    % mic geometry (your ODAS cfgs)
cfg.use_ch  = {1:16, 1:7, [1 15]};                       % 1-based; UMA16B = MIC1 (ch1) + MIC16 (ch15)

% ===================== MEASUREMENTS (replace every NaN) =====================
% For each point: measure X, Y, Z of the mic HOLE with a ruler (cm / 100 = m).
% 'ch' = MATLAB channel of that mic, from the tap test:
%   LAPTOP:  python3 ~/uma16_tap_map.py --device hw:CARD=UMA16A,DEV=0 --mics MIC8 MIC9 MIC2 MIC15
%            python3 ~/uma16_tap_map.py --device <UMA-8 card> --nch 8 --mics <centre label> <ring label a> <ring label b>
%   -> use the printed "MATLAB channel" numbers below.
%
% UMA16A (intact, LEFT): 4 corner mic holes A1..A4
geomA.ch  = [NaN NaN NaN NaN];
geomA.xyz = [ NaN NaN NaN ;      % A1
              NaN NaN NaN ;      % A2
              NaN NaN NaN ;      % A3
              NaN NaN NaN ];     % A4
% UMA8 (top, flat): centre mic + 2 ring mic holes (U8-a, U8-b)
geom8.ch  = [NaN NaN NaN];
geom8.xyz = [ 0   0   NaN ;      % centre = O (only its height Z)
              NaN NaN NaN ;      % U8-a
              NaN NaN NaN ];     % U8-b
% UMA16B (cut, RIGHT): MIC1 and MIC16 holes (B1, B2); channels known: MIC1 = 1, MIC16 = 15
geomB.ch  = [1 15];
geomB.xyz = [ NaN NaN NaN ;      % B1 = MIC1
              NaN NaN NaN ];     % B2 = MIC16
% Sync speaker = where the PHONE lies while it plays sync_chirp.wav (live calibration)
% (offline mode: the speaker that plays the chirp during tri_record)
cfg.spk   = [NaN NaN NaN];
% Side each board listens to: [0 0 1] = mic side UP (MIC labels visible), [0 0 -1] = mic side on the table
cfg.front = {[0 0 1], [0 0 1], [0 0 1]};
% ===========================================================================

cfg.geom  = {geomA, geom8, geomB};
cfg.chirp = struct('f0', 6000, 'f1', 9000, 'dur', 0.03, 'period', 1.0);

% optional: phone at a KNOWN position during [t0 t1] s -> yaw check/correction per board ([] = skip)
cfg.calib   = [];                  % e.g. [5 15  0.0 1.5 0.3]
cfg.autoyaw = true;

% localisation
cfg.seg        = [2 1e9];          % seconds analysed (UMA16A time)
cfg.band       = [200 4000];       % Hz
cfg.nfft       = 2048;
cfg.block      = 0.5;              % s per direction estimate
cfg.grid_az    = -180:4:180;
cfg.grid_el    = 0:6:84;           % sources above the table
cfg.upper_only = true;             % all boards lie flat: keep the upper half-space
cfg.range      = Inf;              % Inf = far source; ~1.5 for a phone ~1.5 m away indoors
cfg.truth      = [];               % optional rows [t0 t1 az el] (deg): az 0 = +Y, +90 = +X
cfg.plot       = true;

% the old pose fields are only used when cfg.geom is empty
cfg.pos = {[0 0 0], [0 0 0], [0 0 0]}; cfg.facing = [0 0 0]; cfg.roll = [0 0 0];
cfg.mount = {'flat', 'flat', 'flat'};

% real time (phone = sync source during calibration, then test sound)
cfg.mode         = 'waves';        % 'waves'        = SEE the waves on every mic, delays, delay-and-sum, phase (live)
                                   % 'tdoa_detect'  = ANY SOUND: noise-learnt detection + direction (freq-mode strategy)
                                   % 'freq_live'    = DRONE by FREQUENCY (harmonic comb) + direction, live
                                   % 'tdoa_live'    = OWN TDOA method in REAL TIME (calibration saved)
                                   % 'tdoa_session' = OWN TDOA method: guided calibration + TEST POINTS + error
                                   % 'live3'   = all 3 boards, triangulation (rough positions below, no chirp)
                                   % 'live16a' = UMA16A alone (no measurements)
                                   % 'live'    = all 3 boards as one array (needs the exact MEASUREMENTS + chirp)
                                   % 'offline' = analyse the tri_record files

% 'live3' triangulation: ROUGH positions (m) of the board CENTRES, a few cm accuracy is enough.
% Same frame as above: O = UMA-8 centre, X right (seen from the laptop), Y toward the UMA-8.
% [x y z]: z = height of the board centre above the UMA-8 mics (needed for the 3-D position)
% Layout of 2026-11-02 (cardboard triangle photo + sketch): UMA-8 at the apex, UMA16A (intact)
% on the LEFT side, UMA16B (cut) on the RIGHT side, laptop at the bottom.  From the sketch:
% UMA-8 -> left board 14 cm, UMA-8 -> right board 15 cm, between the board bottoms 18.5 cm.
% Heights unknown: fill in z (board mic centre height minus UMA-8 mic height) if known.
cfg.tri.pos      = {[-0.073 -0.239  0.00], ...  % UMA16A (LEFT, intact)  <-- check with a ruler
                    [ 0.000  0.000  0.00], ...  % UMA8 (origin, apex)
                    [ 0.076 -0.249  0.00]};     % UMA16B (RIGHT, cut)    <-- check with a ruler
% 2 calibration points where you hold the phone (white noise) at start, ~1 m from the boards.
% Mark them on the floor/table; a few cm accuracy is enough.
cfg.tri.cal      = [ 0.00  1.00 ;         % cal 1: 1 m straight ahead, beyond the UMA-8
                     1.00 -0.20 ];        % cal 2: 1 m to the right
cfg.tri.cal_sec  = 10;                    % s of noise at each point
cfg.tri.settle_sec = 2;                   % first s of each point ignored (phone settling)
cfg.tri.move_sec = 6;                     % s to move the phone before each point
% UMA16B (cut): only MIC1 (ch1) and MIC16 (ch15) remain -> positions of MIC1 and MIC16, 42 mm apart
cfg.tri.b_mics   = [1 16];
cfg.tri.b_weight = 0.5;                   % weight of UMA16B in the position fit
cfg.tri.sigma_deg = 5;                    % expected bearing error of UMA16A / UMA8
cfg.tri.sigma_us  = 10;                   % expected time-difference error of UMA16B
% own TDOA method (tdoa_estimate / tdoa_solve / tdoa_calibrate)
cfg.tdoa.nfft    = 2048;                  % frame for the cross-spectra
cfg.tdoa.band    = [300 6000];            % Hz
cfg.tdoa.bands   = [300 700; 700 1500; 1500 3000; 3000 6000];   % sub-band consensus
cfg.tdoa.up      = 8;                     % lag grid = 1/(8 fs) = 2.6 us
cfg.tdoa.src_z   = 0;                     % source height relative to the mics (m)
cfg.tdoa.sigma   = 10e-6;                 % expected TDOA error (s)
cfg.tdoa.huber   = 2;                     % robust threshold (sigmas)
cfg.tdoa.iters   = 15;                    % Gauss-Newton iterations
cfg.tdoa.board_weight = [1 1 0.5];        % UMA16A, UMA8, UMA16B
cfg.tdoa.adaptive = true;                 % per block, a board whose TDOAs do not fit one direction
                                          % gets less weight: (best residual / its residual)^2
% Calibration points [x y z] (m), z = phone height above the UMA-8 mics.
% Use different heights: they fix the tilt of every board, needed for a good z.
cfg.tdoa.cal     = [ 0.00  1.00  0.00 ;   % cal 1: 1 m ahead, at mic height (on the table)
                     1.00 -0.20  0.00 ;   % cal 2: 1 m to the right, on the table
                    -1.00 -0.20  0.00 ;   % cal 3: 1 m to the left, on the table
                     0.00  1.00  0.40 ];  % cal 4: 1 m ahead, held 40 cm ABOVE the table
%                    1.00 -0.20  0.60 ];  % (cal 5, removed: it broke the 2026-11-02 10:07 calibration)
% Test points [x y z] (m) with KNOWN positions -> tdoa_session measures the error at each.
% z = height above the UMA-8 mics: hold the phone at that height (0 = same height as the boards).
cfg.tdoa.test    = [ 0.00  0.50  0.00 ;   % 1: ahead, 0.5 m
                     0.00  2.00  0.00 ;   % 2: ahead, 2 m
                     0.00  3.00  0.00 ;   % 3: ahead, 3 m
                     0.71  0.71  0.00 ;   % 4: 45 deg right, 1 m
                     1.41  1.41  0.00 ;   % 5: 45 deg right, 2 m
                    -0.71  0.71  0.00 ;   % 6: 45 deg left, 1 m
                    -1.41  1.41  0.00 ;   % 7: 45 deg left, 2 m
                     0.00  1.50  0.50 ;   % 8: ahead, 1.5 m, 50 cm UP
                     0.71  0.71  0.40 ];  % 9: 45 deg right, 1 m, 40 cm UP
cfg.tdoa.flat    = [false false false];   % true for a board known to lie FLAT (UMA16A, UMA8, UMA16B):
                                          % only its rotation is calibrated (the UMA-8 is on a sloped
                                          % wedge, so it is NOT flat here)
% Known mounting (photo of 2026-11-02): which way each board's MIC FACE points, as a direction
% in this frame (az 0 = +Y ahead, +90 = +X right).  The UMA-16s stand along the sides of the
% cardboard triangle (17 deg from the Y axis) and face outward.  The calibration keeps each
% normal within tol degrees of it (front/back sign is free) and fixes the exact UMA-16 scale.
cfg.tdoa.prior = struct( ...
    'normal', {[-cosd(17) sind(17) 0], [0 1 0], [cosd(17) sind(17) 0]}, ...  % UMA16A left, UMA8 front, UMA16B right
    'tol',    {25, 25, 25}, ...
    'scale',  {1, NaN, 1});                   % UMA-8: NaN = fitted until its ring radius is measured
cfg.tdoa.solve_z = true;                  % true = 3-D position (x, y, z); false = x, y at src_z
cfg.tdoa.grid_z  = 0.0:0.02:1.2;          % heights searched (m above the UMA-8 mics)
cfg.tdoa.coarse  = 0.06;                  % 3-D start grid step (m)
cfg.tdoa.max_range = 10;                  % search out to 10 m from the boards ([] = the x-y-z box above)
cfg.tdoa.el_grid = -10:5:80;              % elevations searched (deg above the boards)
cfg.tdoa.move_sec = 8;                    % s to move the phone to the next point
cfg.tdoa.rec_sec  = 10;                   % s recorded at each point
cfg.tdoa.settle_sec = 2;                  % first s of each point ignored
cfg.tdoa.reuse_cal = false;               % tdoa_live: reuse tdoa_calibration.mat (false = recalibrate)
cfg.tdoa.max_rms_us = 25;                 % tdoa_live: fit worse than this -> 'uncertain'
cfg.tdoa.max_std  = 0.25;                 % tdoa_live: uncertainty larger than this (m) -> 'uncertain'
cfg.tri.grid_x   = -1.5:0.02:1.5;         % area searched for the source (m)
cfg.tri.grid_y   = -1.0:0.02:1.5;
cfg.live_cal     = 20;             % s of phone chirps at cfg.spk before localisation starts
cfg.mic_gain     = 10;             % software amplifier on all mic signals (x10 = +20 dB)
cfg.live_gate    = -50;            % dBFS AFTER the amplifier: quieter blocks are not localised
                                   % (= -70 dBFS raw; the quiet room is ~-80 raw = -60 amplified)
                                   % NOT used by 'tdoa_detect': an absolute level gate drops every
                                   % distant sound; tdoa_detect gates on contrast above the learnt noise.
                                   % (mic_gain scales signal AND noise: it adds no SNR.)
cfg.live_seconds = Inf;            % stop after this many s (or close the figure)

% frequency mode (freq_live): find the drone's harmonic comb, then its direction from those bins only
cfg.freq.calibrate  = false;              % true = calibrate the boards first (white noise at cfg.tdoa.cal);
                                          % false = reuse tdoa_calibration.mat (refused if cfg.tri.pos changed)
cfg.freq.f0_range   = [100 3000];         % Hz, fundamental (blade-pass) search range
cfg.freq.f0_step_rel = 0.0005;            % f0 grid step (relative) - fine enough for high harmonics
cfg.freq.max_rate   = 0.02;               % f0 may change by 2 % per frame (85 ms): motor speed changes
cfg.freq.fmax       = 9000;               % Hz, highest harmonic used
cfg.freq.n_harm     = 60;                 % harmonics checked (those <= fmax)
cfg.freq.peak_db    = 4;                  % a harmonic counts when, along its track, it is this far above the noise
cfg.freq.min_contrast = 6;                % ... and NARROW: this many dB above the bins 18-24 Hz away
                                          % (rejects room / box / speaker resonances)
cfg.freq.min_coh_h  = 0.5;                % ... and STEADY: its frame-to-frame phase advance changes smoothly
cfg.freq.min_coh    = 0.7;                % mean steadiness of the harmonics that count (resonance noise ~0.4)
cfg.freq.min_harm   = 3;                  % harmonics needed for a detection (2 of them among h = 1..3)
cfg.freq.cfar_k     = 0;                  % path score threshold: median + cfar_k x spread of all candidates
                                          % 0: a strong comb fills all candidates with near-copies of itself and
                                          % masks itself at 6; min_score + harmonic tests + confirm still apply
cfg.freq.min_score  = 30;                 % ... and at least this (sum of dB of the harmonics found)
cfg.freq.n_cand     = 30;                 % best f0 paths checked per update
cfg.freq.nfft       = 8192;               % 5.9 Hz resolution at 48 kHz
cfg.freq.listen_sec = 2;                  % audio used per update (s)
cfg.freq.hop_sec    = 0.5;                % one update every hop_sec
cfg.freq.confirm    = [3 5];              % DRONE when found in 3 of the last 5 updates ...
cfg.freq.max_f0_jump = 0.1;               % ... with f0 changing by at most 10 % between them
cfg.freq.dir_fmin   = 800;                % Hz: only harmonics above this give the direction
cfg.freq.el_range   = [-10 80];           % elevations searched (deg)
cfg.freq.az_step    = 4;                  % direction grid step (deg), refined to 0.5 deg

% sound detection mode (tdoa_detect): same strategy as freq_live, for any sound
%   1. learn each board's noise auto/cross-spectra (fans, PC, electrical tones)
%   2. subtract it; only bins whose coherence stands above the noise vote (tdoa_estimate)
%   3. all 142 pairs vote together for each direction (tdoa_srp), no per-pair decision first
%   4. detection = vote contrast above a CFAR threshold learnt on noise-only blocks
%   5. confirmed when m of the last n updates agree in direction
cfg.det.source      = 'live';             % 'live' = the 3 boards; or a folder with uma16a/uma8/uma16b.wav
cfg.det.cal_file    = 'tdoa_calibration.mat';   % board poses (from tdoa_session / freq calibration)
cfg.det.block       = 1.0;                % s of audio per estimate (longer = weaker sounds, slower reaction)
cfg.det.hop         = 0.25;               % one update every hop s
cfg.det.learn_sec   = 8;                  % first s: noise learnt, NO test sound during them
cfg.det.noise_alpha = 0.03;               % noise update per noise-only update (~8 s memory at hop 0.25)
cfg.det.cfar_k      = 4;                  % threshold = median + cfar_k x spread of the noise-only contrast
cfg.det.min_z       = 3.5;                % ... and at least this (noise-only contrast ~1.6, p95 ~2.5 in simulation)
cfg.det.hist_len    = 240;                % noise-only updates kept for the threshold (60 s)
cfg.det.confirm     = [3 5];              % SOUND (confirmed) when found in 3 of the last 5 updates ...
cfg.det.max_jump    = 15;                 % ... within this many deg of the current direction
cfg.det.az_step     = 2;                  % direction grid (deg), refined to 0.5 deg
cfg.det.ranges      = [0.5 1 2 5 Inf];    % ranges scanned (m); only < ~1 m is resolved by this aperture
cfg.det.plot        = true;
cfg.det.log_file    = 'tdoa_detect_log.mat';

% wave viewer (cfg.mode = 'waves'): raw waves per mic, wavefront, delay-and-sum, delays, phase
cfg.wave.source  = 'live';                % 'live', or a saved session: 'live_sessions\20261103_072819_*'
cfg.wave.cal_file = 'tdoa_calibration.mat';   % for the predicted delays (red circles); optional
cfg.wave.band    = [300 4000];            % Hz shown (comb test: try [300 2000] for cleaner waves)
cfg.wave.block   = 0.2;                   % s of audio used per update for the delays
cfg.wave.hop     = 0.1;                   % s between updates
cfg.wave.win_ms  = 6;                     % ms of waveform shown (delays are < 0.52 ms)
cfg.wave.nfft    = 1024;                  % frame for the delay estimate
cfg.wave.up      = 8;                     % display oversampling (smooth sub-sample shifts)

if strcmp(cfg.mode, 'waves')
    tri_waves_live(cfg);
elseif strcmp(cfg.mode, 'tdoa_detect')
    lg = tdoa_detect_live(cfg);
elseif strcmp(cfg.mode, 'freq_live')
    freq_live(cfg);
elseif strcmp(cfg.mode, 'tdoa_live')
    tdoa_live(cfg);
elseif strcmp(cfg.mode, 'tdoa_session')
    tdoa_session(cfg);
elseif strcmp(cfg.mode, 'live3')
    tri_live3(cfg);
elseif startsWith(cfg.mode, 'live')
    tri_live(cfg);
else
    res = tri_tdoa_combine(cfg);
end
