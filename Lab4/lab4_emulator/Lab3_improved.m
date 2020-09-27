clear all;
close all;

%% Part 1 - Signal Generator and Acquisition

voltageRange = 1;                   % Accepts ±0.1 V, ±0.2 V, ±0.5 V, ±1 V, ±2 V, ±5 V, ±10 V
samplingFrequency = 9*1000;         % Any frequency from 100 to 1 MHz
samples = 10000;                     % Any length up to 1000000 points
waveform = 'sine';                 % Also accepts 'triangle' and 'square'
frequency = 1*1000;                   % Any frequency up to 10 MHz
amplitude = 1;                      % Any amplitude from 0.1 V to 20 V
phase = 0;                         % Any phase shift from 0 to 2pi radians

%y = signal_generator(voltageRange, samplingFrequency, samples, waveform, frequency, amplitude, phase);


figure(1);
t = (1:samples)/samplingFrequency;
y = 1 + 2*sin(2*pi*10000*t);
subplot(2, 1, 1);
plot(t, y);
tmp = axis;
tmp(3) = -voltageRange;
tmp(4) = voltageRange;
axis(tmp);

Y = fft(y);                  % My function DFT(y) is not suited for such long signals...
DSP = abs(Y)/samples;        % Double Side Spectrum
SSP = DSP(1:samples/2+1);    % Single Side Spectrum
SSP(2:end) = 2*SSP(2:end);   % Renormalize Power

deltaf = samplingFrequency / samples;
f = (0:samples/2) * deltaf;

subplot(2,1,2)
plot(f, SSP)

%% Part 2 - Vibration experiment simulator
voltageRange = 2;                  % Accepts ±0.1 V, ±0.2 V, ±0.5 V, ±1 V, ±2 V, ±5 V, ±10 V
samplingFrequency = 1000;           % Sampling frequency from 100 to 1 MHz
samples = 10000;                     % Any length up to 1000000 points

voltage = 12;                       % Motor power supply from 1 V to 15 V
unbalancing = 100;                  % Motor unbalancing weight from 0 g to 100 g
health = 0;                       % Motor state of health from 0% to 100%

y = vibration_generator(voltageRange, samplingFrequency, samples, voltage, unbalancing, health);
fprintf('Mean of the signal = %.3f\n',mean(y))
fprintf('RMS of the signal = %.3f\n',rms(y))
fprintf('RMS of the signal - mean = %.3f\n',rms(y - mean(y)))

t = (1:samples)/samplingFrequency;
subplot(2, 1, 1);
plot(t, y)
tmp = axis;
tmp(3) = -voltageRange;
tmp(4) = voltageRange;
axis(tmp);

Y = fft(y);                  % My function DFT(y) is not suited for such long signals...
DSP = abs(Y)/samples;        % Double Side Spectrum
SSP = DSP(1:samples/2+1);    % Single Side Spectrum
SSP(2:end) = 2*SSP(2:end);   % Renormalize Power

deltaf = samplingFrequency / samples;
f = (0:samples/2) * deltaf;

subplot(2,1,2)
%plot(f, SSP)
plot(f, 20*log10(SSP))
%% Part 5-7 Data Acquisition
%% Voltage
max_vibration = [];
at_freq = [];
signals = [];
spectrums = [];

 for v = 1:12
     voltage = v;                    % Motor power supply from 1 V to 15 V
     unbalancing = 0;                  % Motor unbalancing weight from 0 g to 100 g
     health = 100;                       % Motor state of health from 0% to 100%

     y = vibration_generator(voltageRange, samplingFrequency, samples, voltage, unbalancing, health);
     signals = [signals; y];
     Y = fft(y);                  % My function DFT(y) is not suited for such long signals...
     DSP = abs(Y)/samples;        % Double Side Spectrum
     SSP = DSP(1:samples/2+1);    % Single Side Spectrum
     SSP(2:end) = 2*SSP(2:end);   % Renormalize Power
     spectrums = [spectrums; SSP];
     
     [val, idx] = findpeaks(20*log10(SSP),'MinPeakHeight',-50); % Look for peak avoiding the continuos component
     max_vibration = [max_vibration val(end)];
     at_freq = [at_freq f(idx(end))];
 end
v = 1:12;
save("voltages.mat",'signals', 'spectrums', 'max_vibration', 'at_freq','v')
figure(3)
%subplot(2,1,1)
%plot(v,max_vibration)
%subplot(2,1,2)
plot(v,at_freq,'o-')

p_v = polyfit(v,at_freq,1);

speed_f = ((at_freq-p_v(2))/p_v(1)) / 50;
speed = at_freq / 50;
%% Weight
max_vibration = [];
at_freq = [];
signals = [];
spectrums = [];

 %for w = 0:10:100
 for v = 1:12
     
     voltage = v;                      % Motor power supply from 1 V to 15 V
     unbalancing = 100;                  % Motor unbalancing weight from 0 g to 100 g
     health = 100;                       % Motor state of health from 0% to 100%

     y = vibration_generator(voltageRange, samplingFrequency, samples, voltage, unbalancing, health);
     signals = [signals; y];
     Y = fft(y);                  % My function DFT(y) is not suited for such long signals...
     DSP = abs(Y)/samples;        % Double Side Spectrum
     SSP = DSP(1:samples/2+1);    % Single Side Spectrum
     SSP(2:end) = 2*SSP(2:end);   % Renormalize Power
     spectrums = [spectrums; SSP];
     
     [val, idx] = max(SSP(f<=100 & f>1)); %findpeaks(20*log10(SSP),'MinPeakHeight',-50,'NPeaks',1);%max(SSP(f<=100 & f>1)); % Look for peak at <10kHz
     max_vibration = [max_vibration val];
     at_freq = [at_freq f(idx)];
 end
w = 0:10:100;
v = 1:12;
save("weight.mat",'signals', 'spectrums', 'max_vibration', 'at_freq','w')
plot(v,max_vibration)
%figure(4)
%subplot(2,1,1)
%plot(w,max_vibration)
%hold on
%p_w = polyfit(w(2:end),max_vibration(2:end),1);
%plot(w(2:end),p_w(1)*w(2:end)+p_w(2))
%hold off
%subplot(2,1,2)
%plot(w,at_freq)
%% Health
max_vibration = [];
at_freq = [];
signals = [];
spectrums = [];

%for v = 1:12
for h = 0:5:100
  
    voltage = 12;                      % Motor power supply from 1 V to 15 V
    unbalancing = 0;                  % Motor unbalancing weight from 0 g to 100 g
    health = h;                     % Motor state of health from 0% to 100%
    
    y = vibration_generator(voltageRange, samplingFrequency, samples, voltage, unbalancing, health);
    signals = [signals; y];
    Y = fft(y);                  % My function DFT(y) is not suited for such long signals...
    DSP = abs(Y)/samples;        % Double Side Spectrum
    SSP = DSP(1:samples/2+1);    % Single Side Spectrum
    SSP(2:end) = 2*SSP(2:end);   % Renormalize Power
    spectrums = [spectrums; SSP];
    
    %[max_value, max_freq] = max(SSP(10:end)); % Main vibration
    %[val, idx] = maxk(SSP(10:max_freq-100),4); % Look for 4 health's peak avoiding the constant and main component
    [upper_bound, upper_freq] = find_max_peak(SSP);
    
    val = [];
    idx = [];
    for i = 1:4
        l_range = round(upper_freq * 0.75 - 25);
        h_range = round(upper_freq * 0.75 + 25);
        [p_value,p_freq]  =  max(SSP(l_range:h_range));
        val = [val; p_value];
        idx = [idx; p_freq + l_range];
        upper_freq = upper_freq * 0.75;
    end
    %[max_value, max_freq] = findpeaks(20*log10(SSP(1:upper_freq-100)),'MinPeakHeight',-61, 'MinPeakDistance',50)
    %[idx, i] = sort(idx);        % Make sure they are in ascending order
    %val = val(i);
    max_vibration = [max_vibration; val.'];
    at_freq = [at_freq; f(idx.')];
end
h = 0:5:100;
%save("health.mat",'signals', 'spectrums', 'max_vibration', 'at_freq','h')
y = mean(max_vibration,2).'
%plot(h,max_vibration)
%hold on
plot(h,mean(max_vibration,2))
%%
figure(5)
subplot(2,1,1)
plot(h,max_vibration)
subplot(2,1,2)
plot(h,at_freq)
mean_vib = mean(max_vibration(1:find(h==70),:),2); % take average for health <= 75%
p_h = polyfit(h(1:find(h==70)),mean_vib.',1);
%% Part 8 - Spectrums Plots
load('voltages.mat','spectrums','v')
figure('Name','Spectrum for different Voltages')
n = size(spectrums,1);
for i = 1:n
    subplot(6,2,i)
    plot(f,20*log10(spectrums(i,:)))
    title(['V=',num2str(v(i)),'V'])
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

load('weight.mat','spectrums','w')
figure('Name','Spectrum for different Weights')
n = size(spectrums,1);
for i = 1:n
    subplot(6,2,i)
    %plot(f,spectrums(i,:))
    plot(f,20*log10(spectrums(i,:)))
    title(['W=',num2str(w(i)),'g'])
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

load('health.mat','spectrums','h')
figure('Name','Spectrum for different Healths')
n = size(spectrums,1);
for i = 1:n
    subplot(7,3,i)
    plot(f,20*log10(spectrums(i,:)))
    title(['H=',num2str(h(i)),'%'])
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%% Part 10 - Features extractions
spectrum = SSP; % Load any spectrum here
figure
plot(f,20*log10(spectrum))

ex_peak = f(find(spectrum==max(spectrum(50:end)))); % Extimated main peak
ex_v = (ex_peak-p_v(2))/p_v(1);
fprintf('Estimated power = %.3f V\n',ex_v)

ex_w_pk = max(spectrum(2:40)); % extimated weight peak at low freq.
ex_w_f = f(find(spectrum==ex_w_pk)); % Extimated weight induced frequency
if ex_w_pk > 3*rms(spectrum(2:end)) % Check if it is statistically significant
    ex_w = (ex_w_pk-p_w(2))/p_w(1);
fprintf('Estimated weight = %.3f g\n',ex_w)
else
    fprintf("Weight too light")
end

range = find(f > (ex_w_f+10) & f < (ex_peak-80) ); % range of frequency to estimate health
ex_h_pk = mean(maxk(spectrum(range),4)); % take the average of 4 top peaks
if ex_h_pk > 3*rms(spectrum(2:end)) % Check if it is statistically significant
    ex_h = (ex_h_pk-p_h(2))/p_h(1);
fprintf('Estimated health = %.3f %\n',ex_h)
else
    fprintf("Health > 75%")
end
%% PREPARATION SCRIPT 
% Here I have tried to understand how to compute the DFT and plot the
% Single Side Power Spectrum

fs = 100;               % sampling frequency
t = (0:(1/fs):(10-1/fs))'; % time vector
S = cos(2*pi*15*t) + 0.3*sin(2*pi*25*t);
N = length(S);
X = DFT(S);
%f = (0:N-1)*(fs/N);     % frequency range
X = fftshift(X);         % center signal
fshift = (-N/2:N/2-1)*(fs/N); % centered frequency range
power = abs(X).^2/N;     % zero-centered power
plot(fshift,power)
axis([0 fs/2 0 inf])
xlabel('Fequency [Hz]')
ylabel('Power')
title('DFT')

% Meand and RMS
x_bar = round(mean(S),4)
x_rms = rms(S)

% DFT
function [Sk] = DFT(sn)
N = length(sn);
n = 0:1:N-1; % row vector for n
k = 0:1:N-1; % row vecor for k
coeff = exp(-1j*2*pi/N); % FT coefficients
nk = n'*k; % creates a NxN matrix of nk values
weights = coeff .^ nk; % DFT matrix
Sk = weights*sn;
end

function SSP = gensignal(voltageRange, samplingFrequency, samples, voltage, unbalancing, health)
y = vibration_generator(voltageRange, samplingFrequency, samples, voltage, unbalancing, health);
Y = fft(y);                  % My function DFT(y) is not suited for such long signals...
DSP = abs(Y)/samples;        % Double Side Spectrum
SSP = DSP(1:samples/2+1);    % Single Side Spectrum
SSP(2:end) = 2*SSP(2:end);   % Renormalize Power
end

function [max_val,max_idx] = find_max_peak(SSP)
[val, idx] = findpeaks(20*log10(SSP),'MinPeakHeight',-50);
max_val = val(end);
max_idx = idx(end);
end
