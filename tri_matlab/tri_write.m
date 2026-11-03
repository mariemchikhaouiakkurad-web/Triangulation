function tri_write(w, x, lost)
% TRI_WRITE  Write one captured frame x, preceded by the samples the driver dropped.
% dsp.AudioFileWriter needs a constant frame size, so the dropped samples are written as
% whole zero frames (drivers drop whole frames), keeping the file on the same timeline.
for k = 1:floor(double(lost) / size(x, 1))
    w(zeros(size(x)));
end
w(x);
end
