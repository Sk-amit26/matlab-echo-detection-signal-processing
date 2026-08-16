clc; clear; close all;
fs = 1000;                
t = 0:1/fs:1;               
f_tx = 50;                 
f_main = 120;               
delay = 0.2;                
attenuation = 0.7;          
noise_power = 0.5;         
N = length(t);             
freq = (0:N-1)*(fs/N);      
% ================= SIGNAL GENERATION =================
tx = sin(2*pi*f_tx*t);                           
main_signal = sin(2*pi*f_main*t);               
echo_signal = [zeros(1,round(delay*fs)) ...
               tx(1:end-round(delay*fs))];     
noise = sqrt(noise_power)*randn(size(t));       
rx = main_signal + attenuation*echo_signal + noise; 

figure;
subplot(5,2,1);
plot(t, tx, 'LineWidth',1.2); grid on;
title('Transmitted Signal (Time Domain)');
xlabel('Time (s)'); ylabel('Amplitude');

subplot(5,2,2);
plot(freq, abs(fft(tx)),'LineWidth',1.2); grid on; xlim([0 fs/2]);
title('Transmitted Signal (Frequency Domain)');
xlabel('Frequency (Hz)'); ylabel('Magnitude');

subplot(5,2,3);
plot(t, main_signal,'LineWidth',1.2); grid on;
title('Main Signal (Time Domain)');
xlabel('Time (s)'); ylabel('Amplitude');

subplot(5,2,4);
plot(freq, abs(fft(main_signal)),'LineWidth',1.2); grid on; xlim([0 fs/2]);
title('Main Signal (Frequency Domain)');
xlabel('Frequency (Hz)'); ylabel('Magnitude');

subplot(5,2,5);
plot(t, rx,'LineWidth',1.2); grid on;
title('Received Signal (Time Domain)');
xlabel('Time (s)'); ylabel('Amplitude');

subplot(5,2,6);
plot(freq, abs(fft(rx)),'LineWidth',1.2); grid on; xlim([0 fs/2]);
title('Received Signal (Frequency Domain)');
xlabel('Frequency (Hz)'); ylabel('Magnitude');

% =================USING CROSS-CORRELATION =================
L = length(rx) + length(main_signal) - 1;
cross_corr = zeros(1,L);

for i = 1:L
    for j = 1:length(rx)
        idx = i - j + 1;
        if idx > 0 && idx <= length(main_signal)
            cross_corr(i) = cross_corr(i) + rx(j) * main_signal(idx);
        end
    end
end

% Peak detection (location of main signal)
[~, max_idx] = max(cross_corr);
detected_start = max(max_idx - length(main_signal) + 1, 1);

copy_len = min(length(main_signal), length(rx) - detected_start + 1);
detected_signal_cc = zeros(1, length(rx));
detected_signal_cc(detected_start:detected_start+copy_len-1) = main_signal(1:copy_len);
DET_FFT_CC = abs(fft(detected_signal_cc));

subplot(5,2,7);
plot(freq, DET_FFT_CC,'LineWidth',1.2); grid on; xlim([0 fs/2]);
title('Detected Signal Spectrum (using Cross-Correlation)');
xlabel('Frequency (Hz)'); ylabel('Magnitude');

% =================  USING DFT =================
X_dft = zeros(1, N);
for k = 0:N-1
    for n = 0:N-1
        X_dft(k+1) = X_dft(k+1) + rx(n+1)*exp(-1j*2*pi*k*n/N);
    end
end
X_mag = abs(X_dft);
subplot(5,2,8);
plot(freq, X_mag, 'LineWidth',1.2); grid on; xlim([0 fs/2]);
title('Detected Main Signal ( using DFT)');
xlabel('Frequency (Hz)'); ylabel('Magnitude');

% =================USING  AUTOCORRELATION =================
L_auto = 2*length(rx) - 1;
auto_corr = zeros(1,L_auto);

for i = 1:L_auto
    for j = 1:length(rx)
        idx = i - j + 1;
        if idx > 0 && idx <= length(rx)
            auto_corr(i) = auto_corr(i) + rx(j)*rx(idx);
        end
    end
end

% Peak detection (autocorrelation)
[~, max_auto] = max(auto_corr);
detected_start_auto = max_auto - length(rx);
detected_signal_auto = zeros(1, length(rx));

if detected_start_auto > 0 && detected_start_auto <= length(rx)
    len_auto = min(length(rx)-detected_start_auto+1, length(main_signal));
    detected_signal_auto(detected_start_auto:detected_start_auto+len_auto-1) = main_signal(1:len_auto);
else
    detected_signal_auto(1:length(main_signal)) = main_signal;
end

DET_FFT_AUTO = abs(fft(detected_signal_auto));

subplot(5,2,9);
plot(freq,DET_FFT_AUTO, 'LineWidth',1.2); grid on;
xlim([0 fs/2]);
title('Detected Main Signal ( using Autocorrelation)');
xlabel('Frequency (Hz)'); ylabel('Magnitude');

