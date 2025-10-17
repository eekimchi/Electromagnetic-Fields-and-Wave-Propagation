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

beta = 1j.*k0.*sqrt(ee);

Zo = 50;

capacitor = 20e-12;
solutionfrequency = 2.*pi.*(f);

Zl = 0 + (1./(1j.*solutionfrequency.*capacitor));
Zl = Zl.*50;

wavelength = 0.12;

% === Import HFSS Data ===

hfssopenre = readmatrix('Open Stub Real Z.csv');
hfssopenre = hfssopenre';

hfssopenim = readmatrix('Open Stub Imag Z.csv');
hfssopenim = hfssopenim';

hfssshortedre = readmatrix('Shorted Stub Real Z.csv');
hfssshortedre = hfssshortedre';

hfssshortedim = readmatrix('Shorted Stub Imag Z.csv');
hfssshortedim = hfssshortedim';

% Open
l = 0.4*wavelength;

% Zinopen = Zo.*((Zl + (1j.*Zo.*tan(beta.*l)))./(Zo + (1j.*Zl.*tan(beta.*l))));

Zinopen = -1j.*Zo.*cot(beta.*l);

figure(1)
subplot(2,1,1)
title('Real Port Zo Parameters','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 14)
ylabel('Impedance, \Omega','FontSize', 14)
xlim([2 6])
% ylim([-1 1])
hold on
plot(fgraph, hfssopenre(2,:), 'LineWidth', 3, 'Color', 'cyan')
hold on
plot(fgraph, real(Zinopen).*-1, 'LineWidth', 3, 'LineStyle', '--', 'Color', 'magenta')
legend('HFSS', 'MATLAB', 'Location', 'best')

subplot(2,1,2)
title('Imaginary Port Zo Parameters','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 14)
ylabel('Impedance, \Omega','FontSize', 14)
xlim([2 6])
%ylim([-1 1])
hold on
plot(fgraph, hfssopenim(2,:), 'LineWidth', 3, 'Color', 'cyan')
hold on
plot(fgraph, imag(Zinopen), 'LineWidth', 3, 'LineStyle', '--', 'Color', 'magenta')
legend('HFSS', 'MATLAB', 'Location', 'best')



% Shorted
l = (0.5 - 0.15)*wavelength;

Zinshorted = Zo.*((Zl + (1j.*Zo.*tan(beta.*l)))./(Zo + (1j.*Zl.*tan(beta.*l))));

figure(2)
subplot(2,1,1)
title('Real Port Zo Parameters','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 14)
ylabel('Impedance, \Omega','FontSize', 14)
xlim([2 6])
%ylim([-1 1])
hold on
plot(fgraph, hfssshortedre(2,:), 'LineWidth', 3, 'Color', 'cyan')
hold on
plot(fgraph, real(Zinshorted).*-1, 'LineWidth', 3, 'LineStyle', '--', 'Color', 'magenta')
legend('HFSS', 'MATLAB', 'Location', 'best')

subplot(2,1,2)
title('Imaginary Port Zo Parameters','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 14)
ylabel('Impedance, \Omega','FontSize', 14)
xlim([2 6])
%ylim([-1 1])
hold on
plot(fgraph, hfssshortedim(2,:), 'LineWidth', 3, 'Color', 'cyan')
hold on
plot(fgraph, imag(Zinshorted), 'LineWidth', 3, 'LineStyle', '--', 'Color', 'magenta')
legend('HFSS', 'MATLAB', 'Location', 'best')

% === HFSS S Parameters ===

hfssopenS = readmatrix('HFSS Open Stub S Parameters.csv');
hfssopenS = hfssopenS';

hfssshortedS = readmatrix('HFSS Shorted Stub S Parameters.csv');
hfssshortedS = hfssshortedS';

figure(3)
title('HFSS S Parameters','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 14)
ylabel('dB','FontSize', 14)
xlim([2 6])
%ylim([-1 1])
hold on
plot(fgraph, hfssopenS(2,:), 'LineWidth', 3, 'Color', 'cyan')
hold on
plot(fgraph, hfssopenS(3,:), 'LineWidth', 3, 'Color', 'magenta')
legend('S11', 'S12', 'Location', 'best')

figure(4)
title('HFSS S Parameters','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 14)
ylabel('dB','FontSize', 14)
xlim([2 6])
%ylim([-1 1])
hold on
plot(fgraph, hfssshortedS(2,:), 'LineWidth', 3, 'Color', 'cyan')
hold on
plot(fgraph, hfssshortedS(3,:), 'LineWidth', 3, 'Color', 'magenta')
legend('S11', 'S12', 'Location', 'best')


