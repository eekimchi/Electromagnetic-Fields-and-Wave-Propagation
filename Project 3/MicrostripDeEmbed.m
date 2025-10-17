clc
clear
close all


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

W = 4.7e-3;
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

% ABCD parameters for the line itself
Aline = cosh(gamma.*l);
Bline = Zc.*sinh(gamma.*l);
Cline = sinh(gamma.*l)./Zc;
Dline = cosh(gamma.*l);

ABCDline(1,1,:) = Aline;
ABCDline(1,2,:) = Bline;
ABCDline(2,1,:) = Cline;
ABCDline(2,2,:) = Dline;

% Impedance at 5 GHz
capacitor = 20e-12;
solutionfrequency = 2.*pi.*(f);

Z = 100 + (1./(1j.*solutionfrequency.*capacitor));

A = 1.*ones(1, 701);
B = Z.*ones(1, 701);
C = 0.*ones(1, 701);
D = 1.*ones(1, 701);

ABCD(1,1,:) = A;
ABCD(1,2,:) = B;
ABCD(2,1,:) = C;
ABCD(2,2,:) = D;

% ABCDtotal = pagemtimes(pagemtimes(ABCDline, ABCD), ABCDline);
% 
% Atotal = ABCDtotal(1,1,:);
% Atotal = squeeze(Atotal);
% Atotal = Atotal';
% Btotal = ABCDtotal(1,2,:);
% Btotal = squeeze(Btotal);
% Btotal = Btotal';
% Ctotal = ABCDtotal(2,1,:);
% Ctotal = squeeze(Ctotal);
% Ctotal = Ctotal';
% Dtotal = ABCDtotal(2,2,:);
% Dtotal = squeeze(Dtotal);
% Dtotal = Dtotal';

% convert ABCD to S

% S11 = (Atotal + (Btotal./Zc) - (Ctotal.*Zc) - Dtotal)./(Atotal + (Btotal./Zc) + (Ctotal.*Zc) + Dtotal);
% S12 = (2.*(Atotal.*Dtotal - Btotal.*Ctotal))./(Atotal + (Btotal./Zc) + (Ctotal.*Zc) + Dtotal);
% S21 = (2)./(Atotal + (Btotal./Zc) + (Ctotal.*Zc) + Dtotal);
% S22 = (-Atotal + (Btotal./Zc) - (Ctotal.*Zc) - Dtotal)./(Atotal + (Btotal./Zc) + (Ctotal.*Zc) + Dtotal);

S11 = (A + (B./Zc) - (C.*Zc) - D)./(A + (B./Zc) + (C.*Zc) + D);
S12 = (2.*(A.*D - B.*C))./(A + (B./Zc) + (C.*Zc) + D);
S21 = (2)./(A + (B./Zc) + (C.*Zc) + D);
S22 = (-A + (B./Zc) - (C.*Zc) - D)./(A + (B./Zc) + (C.*Zc) + D);


S(1,1,:) = S11;
S(1,2,:) = S12;
S(2,1,:) = S21;
S(2,2,:) = S22;

% S11 = squeeze(S11);
% S12 = squeeze(S12);
% S21 = squeeze(S21);
% S22 = squeeze(S22);

% === Import HFSS Data ===

hfssSdB = readmatrix('Lump RLC S dB Scale.csv');
hfssSdB = hfssSdB';

hfssSreal = readmatrix('Lump RLC S Real.csv');
hfssSreal = hfssSreal';

hfssSimag = readmatrix('Lump RLC S Imag.csv');
hfssSimag = hfssSimag';

% Initial Results

figure(1)
title('S Parameters','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 16)
ylabel('dB','FontSize', 16)
xlim([2 6])
%ylim([-400 0])
hold on
plot([fgraph NaN fgraph], [hfssSdB(2,:) NaN hfssSdB(3,:)], 'LineWidth', 3, 'Color', 'cyan')
hold on
plot([fgraph NaN fgraph], [mag2db(abs(S22)) NaN mag2db(abs(S12))], 'LineWidth', 3, 'LineStyle', '--', 'Color', 'magenta')
legend('HFSS', 'MATLAB', 'Location', 'best')


% De-embed

% hfssS11 = hfssSreal(2,:) + 1j.*hfssSimag(2,:);
hfssS12 = hfssSreal(3,:) + 1j.*hfssSimag(3,:);
hfssS21 = hfssSreal(4,:) + 1j.*hfssSimag(4,:);
% hfssS22 = hfssSreal(5,:) + 1j.*hfssSimag(5,:);

hfssS11 = 0.*ones(1,701);
hfssS22 = hfssS11;

hfssA = ((1 + hfssS11).*(1 - hfssS22) + (hfssS12.*hfssS21))./(2.*hfssS21);
hfssB = Zc.*(((1 + hfssS11).*(1 + hfssS22) - (hfssS12.*hfssS21))./(2.*hfssS21));
hfssC = (1./Zc).*(((1 - hfssS11).*(1 - hfssS22) - (hfssS12.*hfssS21))./(2.*hfssS21));
hfssD = ((1 - hfssS11).*(1 + hfssS22) + (hfssS12.*hfssS21))./(2.*hfssS21);

hfss(1,1,:) = hfssA;
hfss(1,2,:) = hfssB;
hfss(2,1,:) = hfssC;
hfss(2,2,:) = hfssD;

ABCDinv = pageinv(ABCDline);

hfssrlc = pagemtimes(pagemtimes(ABCDinv, hfss), ABCDinv);

hfssArlc = hfssrlc(1,1,:);
hfssArlc = squeeze(hfssArlc);
hfssBrlc = hfssrlc(1,2,:);
hfssBrlc = squeeze(hfssBrlc);
hfssCrlc = hfssrlc(2,1,:);
hfssCrlc = squeeze(hfssCrlc);
hfssDrlc = hfssrlc(2,2,:);
hfssDrlc = squeeze(hfssDrlc);

% Myabe the Z parameters work out better

hfssrealz = readmatrix('Lump RLC Z Real.csv');
hfssrealz = hfssrealz';

hfssimagz = readmatrix('Lump RLC Z Imag.csv');
hfssimagz = hfssimagz';


hfssZ11 = hfssrealz(2,:) + 1j.*hfssimagz(2,:);
hfssZ21 = hfssrealz(3,:) + 1j.*hfssimagz(3,:);

hfssA = hfssZ11./hfssZ21;
hfssB = ((hfssZ11.*hfssZ11) - (hfssZ21.*hfssZ21))./(hfssZ21);
hfssC = 1./hfssZ21;
hfssD = hfssZ11./hfssZ21;

% hfss(1,1,:) = hfssA;
% hfss(1,2,:) = hfssB;
% hfss(2,1,:) = hfssC;
% hfss(2,2,:) = hfssD;
% 
% hfssrlc = pagemtimes(pagemtimes(ABCDinv, hfss), ABCDinv);
% 
% hfssArlc = hfssrlc(1,1,:);
% hfssArlc = squeeze(hfssArlc);
% hfssBrlc = hfssrlc(1,2,:);
% hfssBrlc = squeeze(hfssBrlc);
% hfssCrlc = hfssrlc(2,1,:);
% hfssCrlc = squeeze(hfssCrlc);
% hfssDrlc = hfssrlc(2,2,:);
% hfssDrlc = squeeze(hfssDrlc);


% only the B parameters matter for impedance

figure(2)
subplot(1,2,1)
title('Real B Parameters','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 14)
ylabel('Impedance, \Omega','FontSize', 14)
xlim([2 6])
%ylim([-1 1])
hold on
% plot(fgraph, real(hfssBrlc), 'LineWidth', 3, 'Color', 'cyan')
plot(fgraph, real(hfssB), 'LineWidth', 3, 'Color', 'cyan')
hold on
plot(fgraph, real(B), 'LineWidth', 3, 'LineStyle', '--', 'Color', 'magenta')
legend('HFSS', 'MATLAB', 'Location', 'best')

subplot(1,2,2)
title('Imaginary B Parameters','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 14)
ylabel('Impedance, \Omega','FontSize', 14)
xlim([2 6])
%ylim([-1 1])
hold on
% plot(fgraph, imag(hfssBrlc), 'LineWidth', 3, 'Color', 'cyan')
plot(fgraph, imag(hfssB), 'LineWidth', 3, 'Color', 'cyan')
hold on
plot(fgraph, imag(B), 'LineWidth', 3, 'LineStyle', '--', 'Color', 'magenta')
legend('HFSS', 'MATLAB', 'Location', 'best')


