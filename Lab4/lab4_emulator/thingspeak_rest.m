close all;
clear all;
clc;


% Vibration experiment simulator
voltageRange = 2;                   % Accepts ±0.1 V, ±0.2 V, ±0.5 V, ±1 V, ±2 V, ±5 V, ±10 V
samplingFrequency = 1000;           % Sampling frequency from 100 to 1 MHz
samples = 100000;                   % Any length up to 1000000 points

% Generates motor state evolution for speed, unbalancing and gear health (it repeats after 10 minutes)
t = (1:samples)/samplingFrequency;
y = vibration_generator_variable(voltageRange, samplingFrequency, samples);

figure(1);
subplot(2, 1, 1);
plot(t, y)
tmp = axis;
tmp(3) = -voltageRange;
tmp(4) = voltageRange;
axis(tmp);

subplot(2, 1, 2);
Y = fft(y);
DSP = abs(Y/samples);
SSP = DSP(1:1+samples/2);
SSP(2:end) = 2*SSP(2:end);
df = samplingFrequency/samples;
f = (0:(samples/2))*df;
subplot(2, 1, 2);
plot(f, 20*log10(SSP));

% Extract features
[f1,f2,f3] = sensor_processing(y);

% Publish
% as Text
thingSpeakURL = 'http://api.thingspeak.com/';
thingSpeakWriteURL = [thingSpeakURL 'update'];
writeApiKey = 'TN9YI0LXPW61CYSD';
%response = webwrite(thingSpeakWriteURL,'api_key',writeApiKey,'field1',f1,'field2',f2,'field3',f3)

% as JSON
thingSpeakURL = 'http://api.thingspeak.com/';
thingSpeakWriteURL = [thingSpeakURL 'update.json'];
jdata = jdataencode(struct('api_key',writeApiKey,'field1',f1,'field2',f2,'field3',f3));
%response = webwrite(thingSpeakWriteURL,jdata)

% Size comparison
savejson('',struct('api_key',writeApiKey,'field1',f1,'field2',f2,'field3',f3),'entry.json');
savemsgpack('',struct('api_key',writeApiKey,'field1',f1,'field2',f2,'field3',f3),'entry.msgs');

% Continuous upolad
i = 0;
%webwrite(thingSpeakWriteURL,'api_key',writeApiKey,'status','ACTIVE')
while(i < 30) % run for 15 mins
    disp(i)
    y = vibration_generator_variable(voltageRange, samplingFrequency, samples);
    [f1,f2,f3] = sensor_processing(y);
    webwrite(thingSpeakWriteURL,'api_key',writeApiKey,'field1',f1,'field2',f2,'field3',f3)
    pause(30)
    i = i + 1;
end
%webwrite(thingSpeakWriteURL,'api_key',writeApiKey,'status','DISABLED')