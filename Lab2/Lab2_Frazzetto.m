% BDMS - Lab 2: Measurment Uncertainty 
% Paolo Frazzetto

% 2
% generate synthetic data
N = 100;
ex = 0.1;
ey = 0.1;
[x y] = gendata(N,1,1);

dN = [10 25 50 75 100 150 200 500];
for n = 1:length(dN)
   N = dN(n);
   [xdata ydata] = gendata(N,ex,ey);
   [m,b,u_m(1,n),u_b(1,n)] = linreg(xdata,ydata,ex,ey);
end

dex = [0.001 0.01 0.05 0.1 0.5 1 5 10];
N = 100;
for n = 1:length(dex)
   ex = dex(n);
   [xdata ydata] = gendata(N,ex,ey);
   [m,b,u_m(1,n),u_b(1,n)] = linreg(xdata,ydata,ex,ey);
end

dey = [0.001 0.01 0.05 0.1 0.5 1 5 10];
N = 100;
for n = 1:length(dey)
   ey = dey(n);
   [xdata ydata] = gendata(N,ex,ey);
   [m,b,u_m(1,n),u_b(1,n)] = linreg(xdata,ydata,ex,ey);
end
% 3 
input = importdata('lab2_measurements.xlsx'); % Load file
figure(1)

x = input.data(10:16,2);
y = input.data(10:16,3);
N = length(x);

x_err = (x*0.04 + 8)/sqrt(3); % Luxmeter Accuracy: +- (4% rds + 8 dgs)
y_err = abs(y)*0.001/sqrt(3); % Multimeter Accuracy: +- 0.1%

errorbar(x,y,y_err,y_err,x_err,x_err,'o')
xlabel('Luxmeter Readout [Lux]')
ylabel('Sensor Output [V]')
title('Calibration Curve')
grid on
hold on

[m,b,u_m,u_b] = linreg(x,y,x_err,y_err);
lin_reg = m * x + b;
plot(x,lin_reg)
legend('Measurements',['y = (' num2str(m) '+-' num2str(u_m) ') x + (' num2str(b) '+-' num2str(u_b)], 'Location','southeast')
hold off

%4
lum = (y-b)./m
u_lum = sqrt(y_err.^2+u_b.^2./(N.^2)+((x-mean(x))/(sum(x-mean(x)).^2).^2*u_m.^2))

% Monte Carlo Method 
n = 5; % Sample size
p = .95
n_simulations = round(1000 / (1 - p)) % Number of simulations

M_mc =zeros(1,n_simulations);
B_mc =zeros(1,n_simulations);

for j=1:n_simulations
    [x_mc,y_mc] = gendata(n,1,1);
    [M_mc(j),B_mc(j), u_m_mc,u_b_mc]=linreg(x_mc,y_mc,1,1);
end

figure
histogram(M_mc)
title("Monte Carlo - m")
figure
histogram(B_mc)
title("Monte Carlo - b")

m_mc = mean(M_mc)
b_mc = mean(B_mc)

CL = n_simulations*.025; % 95% Confidence Interval
B_mc = sort(B_mc);
M_mc = sort(M_mc);

% Intervals
m_min = m_mc - M_mc(CL)
m_max = M_mc(n_simulations - CL) - m_mc
b_min = b_mc - B_mc(CL)
b_max = B_mc(n_simulations - CL) - b_mc
% These results are reasonable for 95% Confidence Intervals for m = 1.9988
% and b = 1.0020

% Function to compute linear regression given x,y and their uncertainties
function [m,b,u_m,u_b] = linreg(x,y,x_err,y_err)
    N = length(x);
    % Compute linaer regression coefficients and uncertainty
    m = sum( (x-mean(x)).*(y-mean(y)) ) / sum( (x-mean(x)).^2 );
    b = mean(y) - m * mean(x);

    dm_dx = ( (y-mean(y)).*sum((x-mean(x)).^2) - sum( (x-mean(x)).*(y-mean(y)) ).* (2.*(x-mean(x))) ) / (sum( (x-mean(x)).^2 )).^2;
    dm_dy = ( x-mean(x) ) / sum((x-mean(x)).^2);
    u_m = sqrt( sum((dm_dx.^2).*(x_err.^2)) + sum((dm_dy.^2).*(y_err.^2)) );

    db_dx = -mean(x).*(dm_dx) - m/N;
    db_dy = 1/N;
    u_b = sqrt( sum((db_dx.^2).*(x_err.^2)) + sum((db_dy.^2).*(y_err.^2)) );
end

% Function to generate data given N, ex, ey
function [x,y] = gendata(N,ex,ey)
    x = linspace(-5,5,N) + normrnd(0,ex,1,N); % generate synthetic data
    y = 2.*x + 1 + normrnd(0,ey,1,N);
end