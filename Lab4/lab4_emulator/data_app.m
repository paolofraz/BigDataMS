readChID = 1058393; % Features Channel
writeChID = 1065024; % Customer Channel
writeAPIKey = '4P5BXWWEBGX2L03T'; 
   
% Read features
f1 = thingSpeakRead(readChID,'Fields',1);
f2 = thingSpeakRead(readChID,'Fields',2);
f3 = thingSpeakRead(readChID,'Fields',3);
display([f1,f2,f3])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate and publish extrapolated quantities %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Feature 1 - Voltage
% y = p1*x + p2 
% y = feature 1
% x = voltage
%Coefficients:
p1 = 13.889;
p2 = 0.0060606;
ex_v = (f1-p2)/p1; % estimated voltage
speed = ex_v / 50; % Rotation Speed [Hz]
rpm = speed * 60; % Rotation Speed [RPM]

% Feature 2 - Weight
% f2 = (k1*weight)*((k2*speed)^2)
k1 = 0.0008234;
k2 = 4.356; 
if f2 > 0.05 % Check if it is statistically significant above an arbitrary threshold
    weight = round(f2 / (k1*(k2*speed)^2),-1);
else
    weight = 0;
end

% Feature 3 - Health
% f3 = ((k3*(health-100))^2)*(k4*speed);
k3 = 0.0205;
k4 = 0.0085119; % Estimated at 50% health
if f3 > 0.001 % Check if it is statistically significant above an arbitrary threshold
    health = round(100 - sqrt(f3 / (k4*speed*(k3^2))),-1);
else
    health = 100;
end

% Publish
thingSpeakWrite(writeChID,[rpm,weight,health],'WriteKey',writeAPIKey);