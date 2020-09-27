%close all;
clear all;
clc;

% Part 1 - Signal Generator and Acquisition

voltageRange = 2;                   % Accepts ±0.1 V, ±0.2 V, ±0.5 V, ±1 V, ±2 V, ±5 V, ±10 V
samplingFrequency = 100000;         % Any frequency from 100 to 1 MHz
samples = 100000;                     % Any length up to 1000000 points
waveform = 'sine';                % Also accepts 'triangle' and 'square'
frequency = 1000;                   % Any frequency up to 10 MHz
amplitude = 1;                      % Any amplitude from 0.1 V to 20 V
phase = 0;                         % Any phase shift from 0 to 2pi radians

y = signal_generator(voltageRange, samplingFrequency, samples, waveform, frequency, amplitude, phase);

figure(1);
t = (1:samples)/samplingFrequency;
subplot(2, 1, 1);
plot(t, y);
tmp = axis;
tmp(3) = -voltageRange;
tmp(4) = voltageRange;
axis(tmp);

% Part 2 - Vibration experiment simulator
voltageRange = 2;                  % Accepts ±0.1 V, ±0.2 V, ±0.5 V, ±1 V, ±2 V, ±5 V, ±10 V
samplingFrequency = 1000;           % Sampling frequency from 100 to 1 MHz
samples = 100000;                     % Any length up to 1000000 points

voltage = 12;                       % Motor power supply from 1 V to 15 V
unbalancing = 40;                  % Motor unbalancing weight from 0 g to 100 g
health = 50;                       % Motor state of health from 0% to 100%

y = vibration_generator(voltageRange, samplingFrequency, samples, voltage, unbalancing, health);

t = (1:samples)/samplingFrequency;
subplot(2, 1, 2);
plot(t, y)
tmp = axis;
tmp(3) = -voltageRange;
tmp(4) = voltageRange;
axis(tmp);
