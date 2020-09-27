function [f1,f2,f3] = sensor_processing(y)
samplingFrequency = 1000;           % Sampling frequency from 100 to 1 MHz
samples = 100000;                   % Any length up to 1000000 points

Y = fft(y);
deltaf = samplingFrequency / samples;
f = (0:samples/2) * deltaf;
DSP = abs(Y)/samples;        % Double Side Spectrum
SSP = DSP(1:samples/2+1);    % Single Side Spectrum
SSP(2:end) = 2*SSP(2:end);   % Renormalize Power

% Feature 1 - Main peak related to the power (voltage supply) at the
% highest frequency
[val, idx] = max(SSP(find(f==5):end));%findpeaks(20*log10(SSP),'MinPeakHeight',-50); % Look for peak avoiding the continuos and weight components
f1 = f(idx(end)); % Take the highest frequency
upper_freq = idx(end); % Save the index for later use

% Feature 2 - (First)Peak related to the weight at low frequencies
[val, idx] = max(SSP(f<=100 & f>1)); % Look for the first meaningful peak
f2 = val;

% Feature 3 - 4 peaks due to the health condition
[val, idx] = maxk(SSP(find(f==5):upper_freq-50),4); % Search for 4 peaks in the range between the other 2 makor components
% val = [];
% idx = [];
% for i = 1:4 % search peaks over specified intervals
%     l_range = round(upper_freq * 0.75 - 50);
%     h_range = round(upper_freq * 0.75 + 50);
%     [p_value,p_freq]  =  max(SSP(l_range:h_range));
%     val = [val; p_value];
%     idx = [idx; p_freq + l_range];
%     upper_freq = upper_freq * 0.75;
% end
    f3 = mean(val); % take the mean peak amplitude
end