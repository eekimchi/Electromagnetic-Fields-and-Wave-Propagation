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

propagation = readmatrix("Project 2 Propagation.csv");
propagation = propagation';

f = 1.75:0.01:6;

figure(1)
title('Propagation of the Waveguide', 'FontSize', 20)
xlabel('Frequency, GHz','FontSize', 16)
ylabel('Propagation Constant, \beta','FontSize', 16)
xlim([1.75 6])
ylim([0 120])
hold on
plot([f NaN f NaN f NaN f NaN f], [real(beta_zTE10) NaN real(beta_zTE20) NaN real(beta_zTE01) NaN real(beta_zTE11) NaN real(beta_zTE21)], 'LineWidth', 3, 'Color', 'cyan')
hold on
plot([f NaN f NaN f NaN f NaN f], [propagation(2,:) NaN propagation(3,:) NaN propagation(4,:) NaN propagation(5,:) NaN propagation(6,:)], 'LineWidth', 3, 'LineStyle', '--','Color','magenta')
legend('MATLAB','HFSS', 'Location', 'Best')

attenuation = readmatrix("Project 2 Attenuation.csv");
attenuation = attenuation';


figure(2)
title('Attenuation of the Waveguide', 'FontSize', 20)
xlabel('Frequency, GHz','FontSize', 16)
ylabel('Attenuation Constant, \alpha','FontSize', 16)
xlim([1.75 6])
ylim([0 120])
hold on
plot([f NaN f NaN f NaN f NaN f], [imag(beta_zTE10) NaN imag(beta_zTE20) NaN imag(beta_zTE01) NaN imag(beta_zTE11) NaN imag(beta_zTE21)], 'LineWidth', 3, 'Color', 'cyan')
hold on
plot([f NaN f NaN f NaN f NaN f], [attenuation(2,:) NaN attenuation(3,:) NaN attenuation(4,:) NaN attenuation(5,:) NaN attenuation(6,:)], 'LineWidth', 3, 'LineStyle', '--','Color', 'magenta')
legend('MATLAB','HFSS', 'Location', 'Best')


