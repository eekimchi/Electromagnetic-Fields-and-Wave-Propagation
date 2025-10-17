clc
clear
close all

Zc = sqrt(200*50)
% Zc = 200;
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


A = ((Zc)/(60))*sqrt((er + 1)/(2)) + ((er - 1)/(er + 1))*(0.23 + ((0.11)/(er)));

Wd = (8*exp(A))/(exp(2*A) - 2);

W = Wd*1.575;

% Getting  S parameters
Zc = 50;
mu0 = 4*pi*10^(-7);
sigma = 58000000;
er = 2.2;
l = 0.120/50;

f = linspace(2e9, 6e9, 701);
omega = 2.*pi.*f;
fgraph = linspace(2, 6, 701);
c = 3e8;
k0 = (2.*pi.*f)./(c);

W = 4.6e-3;
d = 1.575e-3;

Wd = W/d;

dW = 1/Wd;

ee = ((er + 1)/(2)) + ((er - 1)/(2))*((1)/(sqrt(1 + 12*dW)));

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

% === ABCD Parameters ===

l = 165e-3;

% ABCD parameters for the line itself
A50line = cosh(gamma.*l);
B50line = Zc.*sinh(gamma.*l);
C50line = sinh(gamma.*l)./Zc;
D50line = cosh(gamma.*l);

ABCD50line(1,1,:) = A50line;
ABCD50line(1,2,:) = B50line;
ABCD50line(2,1,:) = C50line;
ABCD50line(2,2,:) = D50line;

W = 1.4116e-3;
d = 1.575e-3;

Wd = W/d;

dW = 1/Wd;

ee = ((er + 1)/(2)) + ((er - 1)/(2))*((1)/(sqrt(1 + 12*dW)));

% The characteristic impedance will be considered as a constant
%   100 ohms for better comparison
Zc = 100;

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

l = 0.12/4;

A100line = cosh(gamma.*l);
B100line = Zc.*sinh(gamma.*l);
C100line = sinh(gamma.*l)./Zc;
D100line = cosh(gamma.*l);

ABCD100line(1,1,:) = A100line;
ABCD100line(1,2,:) = B100line;
ABCD100line(2,1,:) = C100line;
ABCD100line(2,2,:) = D100line;

W = 1.4116e-3;
d = 1.575e-3;

Wd = W/d;

dW = 1/Wd;

ee = ((er + 1)/(2)) + ((er - 1)/(2))*((1)/(sqrt(1 + 12*dW)));

% The characteristic impedance will be considered as a constant
%   100 ohms for better comparison
Zc = 100;

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

l = 165e-3;

A100line = cosh(gamma.*l);
B100line = Zc.*sinh(gamma.*l);
C100line = sinh(gamma.*l)./Zc;
D100line = cosh(gamma.*l);

ABCD100line(1,1,:) = A100line;
ABCD100line(1,2,:) = B100line;
ABCD100line(2,1,:) = C100line;
ABCD100line(2,2,:) = D100line;

W = 0.1674e-3;
d = 1.575e-3;

Wd = W/d;

dW = 1/Wd;

ee = ((er + 1)/(2)) + ((er - 1)/(2))*((1)/(sqrt(1 + 12*dW)));

% The characteristic impedance will be considered as a constant
%   100 ohms for better comparison
Zc = 200;

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

A200line = cosh(gamma.*l);
B200line = Zc.*sinh(gamma.*l);
C200line = sinh(gamma.*l)./Zc;
D200line = cosh(gamma.*l);

ABCD200line(1,1,:) = A200line;
ABCD200line(1,2,:) = B200line;
ABCD200line(2,1,:) = C200line;
ABCD200line(2,2,:) = D200line;

ABCD = pagemtimes(pagemtimes(ABCD50line, ABCD100line), ABCD200line);

A = ABCD(1,1,:);
B = ABCD(1,2,:);
C = ABCD(2,1,:);
D = ABCD(2,2,:);

Zc = 100;

S11 = (A + (B./Zc) - (C.*Zc) - D)./(A + (B./Zc) + (C.*Zc) + D);
S12 = (2.*(A.*D - B.*C))./(A + (B./Zc) + (C.*Zc) + D);
S21 = (2)./(A + (B./Zc) + (C.*Zc) + D);
S22 = (-A + (B./Zc) - (C.*Zc) - D)./(A + (B./Zc) + (C.*Zc) + D);

S11 = squeeze(S11);
S11 = S11';
S12 = squeeze(S12);
S12 = S12';
S21 = squeeze(S21);
S21 = S21';
S22 = squeeze(S22);
S22 = S22';

hfssS = readmatrix('Quarter Wavelength S Parameter.csv');
hfssS = hfssS';


figure(1)
title('S Parameters','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 14)
ylabel('dB','FontSize', 14)
xlim([2 6])
%ylim([-1 1])
hold on
plot(fgraph, hfssS(2,:), 'LineWidth', 3, 'Color', 'cyan')
hold on
plot(fgraph, mag2db(abs(S11)), 'LineWidth', 3, 'LineStyle', '--', 'Color', 'magenta')
legend('HFSS', 'MATLAB', 'Location', 'best')




