function tri_record(seconds)
% TRI_RECORD  Record UMA16A, UMA8 and UMA16B at the same time (Windows: WASAPI/DirectSound)
% while the laptop speaker plays sync_chirp.wav.  Writes the files tri_config.m reads:
%   uma16a.wav, uma8.wav, uma16b.wav
%
%   >> tri_record(60)
%
% Keep the laptop speaker at the position cfg.spk you measured in tri_config.m.
% The chirp sync needs >= 3 chirps (1 per second); 30-60 s is a good length.
if nargin < 1, seconds = 60; end
fs = 48000;
frame = 4800;                                          % 100 ms per read
dev = tri_devices();

if ~isfile('sync_chirp.wav'), sync_chirp_make('sync_chirp.wav', ceil(seconds / 60) + 1); end
[chirp, fsc] = audioread('sync_chirp.wav');
chirp = chirp(1:min(end, round((seconds + 2) * fsc)));

for b = 1:3
    rd{b} = audioDeviceReader('Driver', dev(b).driver, 'Device', dev(b).name, ...
        'NumChannels', dev(b).nch, 'SampleRate', fs, 'SamplesPerFrame', frame, ...
        'OutputDataType', 'single'); %#ok<AGROW>
    setup(rd{b});
    wr{b} = dsp.AudioFileWriter(dev(b).file, 'SampleRate', fs, 'DataType', 'int32'); %#ok<AGROW>
end

% Play on the laptop speaker, not the Windows default output (an MCHStreamer).
outs = audiodevinfo().output;
k = find(contains({outs.Name}, 'Realtek') & contains({outs.Name}, 'Haut-parleurs'), 1);
if isempty(k), error('Laptop speaker (Realtek) not found in audiodevinfo outputs'); end
fprintf('Chirp output: %s\n', outs(k).Name);
player = audioplayer(chirp, fsc, 16, outs(k).ID);
nFrames = ceil(seconds * fs / frame);
overruns = zeros(1, 3);
% Warm-up: the streams start one after another; discard what piled up meanwhile.
for k = 1:10
    for b = 1:3, rd{b}(); end
end
fprintf('Recording %d s ... (chirp playing from the laptop speaker)\n', seconds);
play(player);
for k = 1:nFrames
    for b = 1:3
        [x, nOver] = rd{b}();
        overruns(b) = overruns(b) + nOver;
        if nOver > 0
            fprintf('  %s: %d samples dropped at %.1f s\n', dev(b).file, nOver, k * frame / fs);
        end
        wr{b}(x);
    end
    if mod(k, 10) == 0, fprintf('  %3.0f s\n', k * frame / fs); end
end
stop(player);
for b = 1:3, release(wr{b}); end                       % close the WAV files first
for b = 1:3, release(rd{b}); end
for b = 1:3
    fprintf('%-11s %s', dev(b).file, 'written');
    if overruns(b) > 0
        fprintf('   WARNING: %d samples dropped - close other audio apps and record again', overruns(b));
    end
    fprintf('\n');
end
fprintf('Done. Now run: tri_config\n');
end
