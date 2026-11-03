function dev = tri_devices()
% TRI_DEVICES  Capture devices in the order of cfg.names (UMA16A, UMA8, UMA16B).
%
% Windows renames the two MCHStreamers ("MCHStreamer", "2- MCHStreamer", "4- ...") when
% they are re-plugged, so the boards are identified by their signal instead of their name:
% the intact UMA16A has 16 live channels, the cut UMA16B only a few.
%
% Each board is opened at 48 kHz with WASAPI if possible.  WASAPI refuses when Windows has the
% board set to another rate (e.g. 44.1 kHz after re-plugging); DirectSound is then used instead,
% since it converts the rate inside Windows.
fs = 48000;
r = audioDeviceReader('Driver', 'WASAPI');
wasapi = string(getAudioDevices(r));
r = audioDeviceReader('Driver', 'DirectSound');
dsound = string(getAudioDevices(r));

mch = unique([wasapi(contains(wasapi, 'MCHStreamer')); dsound(contains(dsound, 'MCHStreamer'))], 'stable');
u8 = dsound(contains(dsound, 'miniDSP micArray'));
if numel(mch) ~= 2, error('Expected 2 MCHStreamer (UMA-16) inputs, found %d: plug both in', numel(mch)); end
if isempty(u8), error('UMA-8 (miniDSP micArray) input not found: plug it in'); end

live = zeros(1, 2); drv = strings(1, 2);
for i = 1:2
    [x, drv(i)] = probe(mch(i), 16, fs, wasapi, dsound);
    live(i) = sum(rms(x(4801:end, :)) > 0);
end
[~, a] = max(live);                                   % most live channels = intact UMA16A
b = 3 - a;
[~, drv8] = probe(u8(1), 8, fs, string([]), dsound); % the UMA-8 refuses WASAPI in MATLAB
fprintf('UMA16A = "%s" (%d live ch, %s)   UMA16B = "%s" (%d live ch, %s)   UMA8 = "%s" (%s)\n', ...
    mch(a), live(a), drv(a), mch(b), live(b), drv(b), u8(1), drv8);

dev = struct( ...
    'name',   {char(mch(a)), char(u8(1)), char(mch(b))}, ...
    'nch',    {16, 8, 16}, ...
    'driver', {char(drv(a)), char(drv8), char(drv(b))}, ...
    'file',   {'uma16a.wav', 'uma8.wav', 'uma16b.wav'});
end

function [x, driver] = probe(name, nch, fs, wasapi, dsound)
% open the device at fs, WASAPI first then DirectSound; return 0.6 s of audio and the driver used
for driver = ["WASAPI", "DirectSound"]
    if driver == "WASAPI" && ~any(wasapi == name), continue; end
    if driver == "DirectSound" && ~any(dsound == name), continue; end
    try
        p = audioDeviceReader('Driver', char(driver), 'Device', char(name), 'NumChannels', nch, ...
            'SampleRate', fs, 'SamplesPerFrame', 4800, 'OutputDataType', 'double');
        x = [];
        for k = 1:6, x = [x; p()]; end %#ok<AGROW>
        release(p);
        return
    catch err
        last = err;
    end
end
error('Cannot open "%s" at %d Hz with WASAPI or DirectSound: %s', name, fs, last.message);
end
