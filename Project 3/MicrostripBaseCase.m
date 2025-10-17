clc
clear
close all

Zc = 50;
mu0 = 4*pi*10^(-7);
e0 = 8.85e-12;
sigma = 58000000;
er = 2.2;
l = 0.360;

f = linspace(2e9, 6e9, 701);
omega = 2.*pi.*f;
fgraph = linspace(2, 6, 701);
c = 3e8;
k0 = (2.*pi.*f)./(c);

% ===== this is to initially design the microstrip in hfss ======
% taken from Pozar
B = (377*pi)/(2*Zc*sqrt(er));

Wd = (2/pi)*(B - 1 - log(2*B - 1) + ((er - 1)/(2*er))*(log(B - 1) + 0.39 - ((0.61)/(er))));

% 1.575 is the substrate thickness
W = Wd*1.575;

a = 1/(2*pi);

% this is to check the final hfss dimensions

W = 4.7e-3;
d = 1.575e-3;

Wd = W/d;

dW = 1/Wd;

ee = ((er + 1)/(2)) + ((er - 1)/(2))*((1)/(sqrt(1 + 12*dW)));

Zo = (120*pi)/(sqrt(ee)*(Wd + 1.393 + 0.667*log(Wd + 1.444)));
Zo = Zo.*ones(701,1);


% Balanis
w = 4.7e-3;
h = 1.575e-3;
t = 17.5e-6;

weffh = (w/h) + ((1.25/pi)*(t/h)*(1 + log((2*h)/t)));

hweff = 1/weffh;

ereff = ((er + 1)/2) + (((er-1)/2)*((1 + 12*weffh))^(-1/2));

Zcbalanis = ((120*pi)/sqrt(ereff))/(weffh + 1.393 + 0.667*log(weffh + 1.444));
Zcbalanis = Zcbalanis.*ones(701,1);




% =====================================================

% The characteristic impedance will be considered as a constant
%   50 ohms for better comparison
Zc = 50;

% loss tangent is taken from datasheet
% tand = (er*(ee-1))/(ee*(er-1));
tand = 0.0009;

% attenuation due to dielecric loss
alpha_dielectric = (k0.*er.*(ee-1).*tand)./(2.*sqrt(ee).*(er-1));

Rs = sqrt((omega.*mu0)./(2.*sigma));

% attenuation due to conductor loss
alpha_conductor = Rs./(Zc.*W);

% total attenuation
alpha = alpha_dielectric + alpha_conductor;

% propagation
beta = 1j.*k0.*sqrt(ee);

% gamma
gamma = alpha + beta;

% === Import HFSS Data ===

hfssgamma = readmatrix('50 Ohm Base Case Gamma Plot.csv');
hfssgamma = hfssgamma';

hfssrealz = readmatrix('50 Ohm Base Case Real Z Parameter Plot.csv');
hfssrealz = hfssrealz';

hfssimagz = readmatrix('50 Ohm Base Case Imaginary Z Parameter Plot.csv');
hfssimagz = hfssimagz';

hfssS = readmatrix('50 Ohm Base Case S Parameter Plot.csv');
hfssS = hfssS';

hfssZc = readmatrix('50 Ohm Base Case Port Zo Plot.csv');
hfssZc = hfssZc';

% === Characteristic Impedance ===

figure(1)
title('Characteristic Impedance','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 16)
ylabel('Impedance, \Omega','FontSize', 16)
xlim([2 6])
ylim([45 55])
hold on
plot(fgraph, hfssZc(2,:), 'LineWidth', 3, 'Color', 'cyan')
hold on
plot(fgraph, Zo, 'LineWidth', 3, 'LineStyle', '--', 'Color', 'magenta')
hold on
plot(fgraph, Zcbalanis, 'LineWidth', 3, 'LineStyle', '--', 'Color', 'green')
legend('HFSS', 'MATLAB Pozar','MATLAB Balanis', 'Location', 'best')

% === Gamma ===
figure(2)
subplot(1,2,1)
title('Attenuation','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 16)
ylabel('Attenuation Constant, Np/m','FontSize', 16)
xlim([2 6])
ylim([0 0.3])
hold on
plot(fgraph, hfssgamma(2,:), 'LineWidth', 3, 'Color', 'cyan')
hold on
plot(fgraph, real(gamma), 'LineWidth', 3, 'LineStyle', '--', 'Color', 'magenta')
legend('HFSS', 'MATLAB', 'Location', 'best')

subplot(1,2,2)
title('Propagation','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 16)
ylabel('Propagation Constant, rad/m','FontSize', 16)
xlim([2 6])
ylim([0 180])
hold on
plot(fgraph, hfssgamma(3,:), 'LineWidth', 3, 'Color', 'cyan')
hold on
plot(fgraph, imag(gamma), 'LineWidth', 3, 'LineStyle', '--', 'Color', 'magenta')
legend('HFSS', 'MATLAB', 'Location', 'best')

% === ABCD Parameters ===
A = cosh(gamma.*l);
B = Zc.*sinh(gamma.*l);
C = sinh(gamma.*l)./Zc;
D = cosh(gamma.*l);

ABCD(1,1,:) = A;
ABCD(1,2,:) = B;
ABCD(2,1,:) = C;
ABCD(2,2,:) = D;

hfssZ11 = hfssrealz(2,:) + 1j.*hfssimagz(2,:);
hfssZ21 = hfssrealz(3,:) + 1j.*hfssimagz(3,:);

hfssA = hfssZ11./hfssZ21;
hfssB = ((hfssZ11.*hfssZ11) - (hfssZ21.*hfssZ21))./(hfssZ21);
hfssC = 1./hfssZ21;
hfssD = hfssZ11./hfssZ21;

figure(3)
subplot(3,2,1)
title('Real A Parameters','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 14)
ylabel('Magnitude','FontSize', 14)
xlim([2 6])
ylim([-1 1])
hold on
plot(fgraph, real(hfssA), 'LineWidth', 3, 'Color', 'cyan')
hold on
plot(fgraph, real(A), 'LineWidth', 3, 'LineStyle', '--', 'Color', 'magenta')
legend('HFSS', 'MATLAB', 'Location', 'best')

subplot(3,2,2)
title('Imaginary A Parameters','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 14)
ylabel('Magnitude','FontSize', 14)
xlim([2 6])
ylim([-0.1 0.1])
hold on
plot(fgraph, imag(hfssA), 'LineWidth', 3, 'Color', 'cyan')
hold on
plot(fgraph, imag(A), 'LineWidth', 3, 'LineStyle', '--', 'Color', 'magenta')
legend('HFSS', 'MATLAB', 'Location', 'best')

subplot(3,2,3)
title('Real B Parameters','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 14)
ylabel('Resistance, \Omega','FontSize', 14)
xlim([2 6])
ylim([-20 15])
hold on
plot(fgraph, real(hfssB), 'LineWidth', 3, 'Color', 'cyan')
hold on
plot(fgraph, real(B), 'LineWidth', 3, 'LineStyle', '--', 'Color', 'magenta')
legend('HFSS', 'MATLAB', 'Location', 'best')

subplot(3,2,4)
title('Imaginary B Parameters','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 14)
ylabel('Resistance, \Omega','FontSize', 14)
xlim([2 6])
ylim([-60 60])
hold on
plot(fgraph, imag(hfssB), 'LineWidth', 3, 'Color', 'cyan')
hold on
plot(fgraph, imag(B), 'LineWidth', 3, 'LineStyle', '--', 'Color', 'magenta')
legend('HFSS', 'MATLAB', 'Location', 'best')

subplot(3,2,5)
title('Real C Parameters','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 14)
ylabel('Conductance, 1/\Omega','FontSize', 14)
xlim([2 6])
ylim([-1.5e-3 1.5e-3])
hold on
plot(fgraph, real(hfssC), 'LineWidth', 3, 'Color', 'cyan')
hold on
plot(fgraph, real(C), 'LineWidth', 3, 'LineStyle', '--', 'Color', 'magenta')
legend('HFSS', 'MATLAB', 'Location', 'best')

subplot(3,2,6)
title('Imaginary C Parameters','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 14)
ylabel('Conductance, 1/\Omega','FontSize', 14)
xlim([2 6])
ylim([-0.025 0.025])
hold on
plot(fgraph, imag(hfssC), 'LineWidth', 3, 'Color', 'cyan')
hold on
plot(fgraph, imag(C), 'LineWidth', 3, 'LineStyle', '--', 'Color', 'magenta')
legend('HFSS', 'MATLAB', 'Location', 'best')

% === S Parameters ===
S11 = (A + (B./Zc) - (C.*Zc) - D)./(A + (B./Zc) + (C*Zc) + D);
S12 = (2.*(A.*D - B.*C))./(A + (B./Zc) + (C*Zc) + D);
S21 = (2)./(A + (B./Zc) + (C*Zc) + D);
S22 = (-A + (B./Zc) - (C.*Zc) - D)./(A + (B./Zc) + (C*Zc) + D);

S(1,1,:) = S11;
S(1,2,:) = S12;
S(2,1,:) = S21;
S(2,2,:) = S22;

figure(4)
title('S Parameters','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 16)
ylabel('dB','FontSize', 16)
xlim([2 6])
ylim([-400 0])
hold on
plot([fgraph NaN fgraph], [hfssS(2,:) NaN hfssS(3,:)], 'LineWidth', 3, 'Color', 'cyan')
hold on
plot([fgraph NaN fgraph], [mag2db(abs(S11)) NaN mag2db(abs(S12))], 'LineWidth', 3, 'LineStyle', '--', 'Color', 'magenta')
legend('HFSS', 'MATLAB', 'Location', 'best')



% === Z Parameters ===
Z11 = A./C;
Z12 = (A.*D - B.*C)./C;
Z21 = 1./C;
Z22 = D./C;

Z(1,1,:) = Z11;
Z(1,2,:) = Z12;
Z(2,1,:) = Z21;
Z(2,2,:) = Z22;

figure(5)
subplot(2,2,1)
title('Real Z11 Parameters','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 16)
ylabel('Impedance, \Omega','FontSize', 16)
xlim([2 6])
ylim([0 1600])
hold on
plot(fgraph,hfssrealz(2,:), 'LineWidth', 3, 'Color', 'cyan')
hold on
plot(fgraph,real(Z11), 'LineWidth', 3, 'LineStyle', '--', 'Color', 'magenta')
legend('HFSS', 'MATLAB', 'Location', 'best')

subplot(2,2,2)
title('Imaginary Z11 Parameters','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 16)
ylabel('Impedance, \Omega','FontSize', 16)
xlim([2 6])
ylim([-1000 1000])
hold on
plot(fgraph, hfssimagz(2,:), 'LineWidth', 3, 'Color', 'cyan')
hold on
plot(fgraph, imag(Z11), 'LineWidth', 3, 'LineStyle', '--', 'Color', 'magenta')
legend('HFSS', 'MATLAB', 'Location', 'best')

subplot(2,2,3)
title('Real Z21 Parameters','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 16)
ylabel('Impedance, \Omega','FontSize', 16)
xlim([2 6])
ylim([-2000 2000])
hold on
plot(fgraph,hfssrealz(3,:), 'LineWidth', 3, 'Color', 'cyan')
hold on
plot(fgraph,real(Z21), 'LineWidth', 3, 'LineStyle', '--', 'Color', 'magenta')
legend('HFSS', 'MATLAB', 'Location', 'best')

subplot(2,2,4)
title('Imaginary Z21 Parameters','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 16)
ylabel('Impedance, \Omega','FontSize', 16)
xlim([2 6])
ylim([-1000 1000])
hold on
plot(fgraph, hfssimagz(3,:), 'LineWidth', 3, 'Color', 'cyan')
hold on
plot(fgraph, imag(Z21), 'LineWidth', 3, 'LineStyle', '--', 'Color', 'magenta')
legend('HFSS', 'MATLAB', 'Location', 'best')
