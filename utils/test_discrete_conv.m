clear all;
close all;

signal  = randi([2 8], 16,1);
signalbkl = zeros(24,1);
signalbkl(1:8) = signal(end-7:end);
signalbkl(9:end) = signal;

signalbkl_upsampling = upsample(signalbkl,4);
channel = rand(8,1);

out1 = filter(channel,1,signalbkl_upsampling);

fft_signal = fft(signal);
fft_signal_up = fft(signalbkl_upsampling);
fft_channel = fft(channel,length(signalbkl_upsampling));

out2 = ifft(fft_signal_up.*fft_channel,length(signalbkl_upsampling));
    
plot(out1(9:end));
%hold on
figure
plot(out2);
