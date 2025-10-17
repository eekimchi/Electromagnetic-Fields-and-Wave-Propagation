clc
clear
close all

% ========== WAVEGUIDE ==========

% Waveguide dimensions are for fc = 2 GHz

a = 0.075;
b = 0.0375;

% Yes Dr. Ruyle, this shitty code does work

% Waveguide Modes

% Transverse Electric TE

for m = 1:1:5
for n = 1:1:5
    fcTE(m,n) = ((3e8)./(2.*pi)).*sqrt((((m-1).*pi)./a).^2 + (((n-1).*pi)./b).^2);
end
end

% Transverse Magnetic TM

for m = 1:1:5
for n = 1:1:5
    fcTM(m,n) = ((3e8)./(2.*pi)).*sqrt((((m).*pi)./a).^2 + (((n).*pi)./b).^2);
end
end

% Propagation and Attenuation Constants

f = 1.75e9:0.01e9:6e9;

omega_f = 2.*pi.*f;

omega_fcTE10 = 2*pi*(fcTE(2,1));
omega_fcTE20 = 2*pi*(fcTE(3,1));
omega_fcTE01 = 2*pi*(fcTE(1,2));
omega_fcTE11 = 2*pi*(fcTE(2,2));
omega_fcTE21 = 2*pi*(fcTE(3,2));

beta_zTE10 = sqrt((((omega_f).^2)./(3e8)^2) - (((omega_fcTE10).^2)./(3e8)^2));
beta_zTE20 = sqrt((((omega_f).^2)./(3e8)^2) - (((omega_fcTE20).^2)./(3e8)^2));
beta_zTE01 = sqrt((((omega_f).^2)./(3e8)^2) - (((omega_fcTE01).^2)./(3e8)^2));
beta_zTE11 = sqrt((((omega_f).^2)./(3e8)^2) - (((omega_fcTE11).^2)./(3e8)^2));
beta_zTE21 = sqrt((((omega_f).^2)./(3e8)^2) - (((omega_fcTE21).^2)./(3e8)^2));

% mu
mu = 4*pi*(10^-7);

% the actual port impedances

Z10 = 1j.*(omega_fcTE10.*mu)./beta_zTE10;
Z20 = 1j.*(omega_fcTE20.*mu)./beta_zTE20;
Z01 = 1j.*(omega_fcTE01.*mu)./beta_zTE01;
Z11 = 1j.*(omega_fcTE11.*mu)./beta_zTE11;
Z21 = 1j.*(omega_fcTE21.*mu)./beta_zTE21;

imZo = readmatrix("Project 2 Imaginary Zo.csv");
imZo = imZo';

figure(1)
title('Imaginary Port Impedance','FontSize', 20)
xlabel('Frequency, GHz','FontSize', 16)
ylabel('Impedance','FontSize', 16)
xlim([1.75e9 6e9])
ylim([0 3500])
hold on
plot([f NaN f NaN f NaN f NaN f], [real(Z10) NaN real(Z20) NaN real(Z01) NaN real(Z11) NaN real(Z21)], 'LineWidth', 3, 'Color', 'cyan')
hold on 
plot([f NaN f NaN f NaN f NaN f], [imZo(2,:) NaN imZo(3,:) NaN imZo(4,:) NaN imZo(5,:) NaN imZo(6,:)], 'LineWidth', 3, 'LineStyle', '--','Color','magenta')
legend('MATLAB', 'HFSS', 'Location', 'best')

reZo = readmatrix("Project 2 Real Zo.csv");
reZo = reZo';

figure(2)
title('Real Port Impedance','FontSize', 20)
xlabel('Frequency','FontSize', 16)
ylabel('Impedance','FontSize', 16)
xlim([1.75e9 6e9])
ylim([0 3500])
hold on
plot([f NaN f NaN f NaN f NaN f], [imag(Z10) NaN imag(Z20) NaN imag(Z01) NaN imag(Z11) NaN imag(Z21)], 'LineWidth', 3, 'Color', 'cyan')
hold on 
plot([f NaN f NaN f NaN f NaN f], [reZo(2,:) NaN reZo(3,:) NaN reZo(4,:) NaN reZo(5,:) NaN reZo(6,:)] , 'LineWidth', 3, 'LineStyle', '--','Color','magenta')
legend('MATLAB','HFSS','Location', 'best')
