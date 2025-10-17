clc
clear
close all

% permittivity
e0 = 8.854e-12;

% permeability
mu0 = (4*pi)*(10^(-7));

% frequency range
f = linspace(1.75e9, 6e9, 30);

% convert frequency to omega
omega_f = 2.*pi.*f;

% Waveguide Dimensions
a = 0.075;
b = 0.0375;

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

omega_fcTE10 = 2*pi*(fcTE(2,1));
omega_fcTE20 = 2*pi*(fcTE(3,1));
omega_fcTE01 = 2*pi*(fcTE(1,2));
omega_fcTE11 = 2*pi*(fcTE(2,2));
omega_fcTE21 = 2*pi*(fcTE(3,2));

% beta_z
beta_zTE10 = sqrt((((omega_f).^2)./(3e8).^2) - (((omega_fcTE10).^2)./(3e8).^2));
beta_zTE20 = sqrt((((omega_f).^2)./(3e8).^2) - (((omega_fcTE20).^2)./(3e8).^2));
beta_zTE01 = sqrt((((omega_f).^2)./(3e8).^2) - (((omega_fcTE01).^2)./(3e8).^2));
beta_zTE11 = sqrt((((omega_f).^2)./(3e8).^2) - (((omega_fcTE11).^2)./(3e8).^2));
beta_zTE21 = sqrt((((omega_f).^2)./(3e8).^2) - (((omega_fcTE21).^2)./(3e8).^2));

% beta_x

for m = 1:1:5
    beta_x(m) = ((m-1).*pi)./a;
end

% beta_y

for n = 1:1:5
    beta_y(n) = ((n-1).*pi)./b;
end

% k_c

kc10 = sqrt(beta_x(2).^2 + beta_y(1).^2);
kc20 = sqrt(beta_x(3).^2 + beta_y(1).^2);
kc01 = sqrt(beta_x(1).^2 + beta_y(2).^2);
kc11 = sqrt(beta_x(2).^2 + beta_y(2).^2);
kc21 = sqrt(beta_x(3).^2 + beta_y(2).^2);

% establish coordinates

x = linspace(0, 0.075, 30);
y = linspace(0, 0.0375, 30);
y = y';
z = linspace(0, 0, 30);

% Create a meshgrid for plotting
[X,Y] = meshgrid(x,y);

% E field in X

E_x10 = ((1j.*omega_fcTE10.*mu0.*0.*pi)./(((kc10).^2).*0.0375)).*cos(beta_x(2).*X).*sin(beta_y(1).*Y).*exp(-1j.*beta_zTE10.*z);
E_x20 = ((1j.*omega_fcTE20.*mu0.*0.*pi)./(((kc20).^2).*0.0375)).*cos(beta_x(3).*X).*sin(beta_y(1).*Y).*exp(-1j.*beta_zTE20.*z);
E_x01 = ((1j.*omega_fcTE01.*mu0.*1.*pi)./(((kc01).^2).*0.0375)).*cos(beta_x(1).*X).*sin(beta_y(2).*Y).*exp(-1j.*beta_zTE01.*z);
E_x11 = ((1j.*omega_fcTE11.*mu0.*1.*pi)./(((kc11).^2).*0.0375)).*cos(beta_x(2).*X).*sin(beta_y(2).*Y).*exp(-1j.*beta_zTE11.*z);
E_x21 = ((1j.*omega_fcTE21.*mu0.*1.*pi)./(((kc21).^2).*0.0375)).*cos(beta_x(3).*X).*sin(beta_y(2).*Y).*exp(-1j.*beta_zTE21.*z);

% E field in Y

E_y10 = ((-1j.*omega_fcTE10.*mu0.*1.*pi)./(((kc10).^2).*0.075)).*sin(beta_x(2).*X).*cos(beta_y(1).*Y).*exp(-1j.*beta_zTE10.*z);
E_y20 = ((-1j.*omega_fcTE20.*mu0.*2.*pi)./(((kc20).^2).*0.075)).*sin(beta_x(3).*X).*cos(beta_y(1).*Y).*exp(-1j.*beta_zTE20.*z);
E_y01 = ((-1j.*omega_fcTE01.*mu0.*0.*pi)./(((kc01).^2).*0.075)).*sin(beta_x(1).*X).*cos(beta_y(2).*Y).*exp(-1j.*beta_zTE01.*z);
E_y11 = ((-1j.*omega_fcTE11.*mu0.*1.*pi)./(((kc11).^2).*0.075)).*sin(beta_x(2).*X).*cos(beta_y(2).*Y).*exp(-1j.*beta_zTE11.*z);
E_y21 = ((-1j.*omega_fcTE21.*mu0.*2.*pi)./(((kc21).^2).*0.075)).*sin(beta_x(3).*X).*cos(beta_y(2).*Y).*exp(-1j.*beta_zTE21.*z);

% H Field in X

H_x10 = sin(beta_x(2).*X).*cos(beta_y(1).*Y).*exp(-1j.*beta_zTE10.*z);
H_x20 = sin(beta_x(3).*X).*cos(beta_y(1).*Y).*exp(-1j.*beta_zTE20.*z);
H_x01 = sin(beta_x(1).*X).*cos(beta_y(2).*Y).*exp(-1j.*beta_zTE01.*z);
H_x11 = ((omega_fcTE11.*mu0.*1.*pi)./(((kc11).^2).*0.075)).*sin(beta_x(2).*X).*cos(beta_y(2).*Y).*exp(-1j.*beta_zTE11.*z);
H_x21 = sin(beta_x(3).*X).*cos(beta_y(2).*Y).*exp(-1j.*beta_zTE21.*z);

% H Field in Y

H_y10 = cos(beta_x(2).*X).*sin(beta_y(1).*Y).*exp(-1j.*beta_zTE10.*z);
H_y20 = cos(beta_x(3).*X).*sin(beta_y(1).*Y).*exp(-1j.*beta_zTE20.*z);
H_y01 = cos(beta_x(1).*X).*sin(beta_y(2).*Y).*exp(-1j.*beta_zTE01.*z);
H_y11 = ((omega_fcTE11.*mu0.*1.*pi)./(((kc11).^2).*0.0375)).*cos(beta_x(2).*X).*sin(beta_y(2).*Y).*exp(-1j.*beta_zTE11.*z);
H_y21 = cos(beta_x(3).*X).*sin(beta_y(2).*Y).*exp(-1j.*beta_zTE21.*z);

% Plot E Field

% TE10

figure(1)
quiverC2D(X,Y,imag(E_x10),imag(E_y10))
xlim([0, 0.075])
ylim([0, 0.0375])
title('Waveguide TE_{10} E Field at the Port', 'FontSize', 20)
xlabel('x length, m', 'FontSize', 16)
ylabel('z length, m', 'FontSize', 16)

% TE 20

figure(2)
quiverC2D(X,Y,imag(E_x20),imag(E_y20))
xlim([0, 0.075])
ylim([0, 0.0375])
title('Waveguide TE_{20} E Field at the Port', 'FontSize', 20)
xlabel('x length, m', 'FontSize', 16)
ylabel('z length, m', 'FontSize', 16)

% TE 01

figure(3)
quiverC2D(X,Y,imag(E_x01),imag(E_y01))
xlim([0, 0.075])
ylim([0, 0.0375])
title('Waveguide TE_{01} E Field at the Port', 'FontSize', 20)
xlabel('x length, m', 'FontSize', 16)
ylabel('z length, m', 'FontSize', 16)

% TE 11

figure(4)
quiverC2D(X,Y,imag(E_x11),imag(E_y11))
xlim([0, 0.075])
ylim([0, 0.0375])
title('Waveguide TE_{11} E Field at the Port', 'FontSize', 20)
xlabel('x length, m', 'FontSize', 16)
ylabel('z length, m', 'FontSize', 16)

% TE 21

figure(5)
quiverC2D(X,Y,imag(E_x21),imag(E_y21))
xlim([0, 0.075])
ylim([0, 0.0375])
title('Waveguide TE_{21} E Field at the Port', 'FontSize', 20)
xlabel('x length, m', 'FontSize', 16)
ylabel('z length, m', 'FontSize', 16)

% Plot H Field

% TE 10

figure(6)
quiverC2Dfm(X,Y,real(H_x10),real(H_y10))
xlim([0, 0.075])
ylim([0, 0.0375])
title('Waveguide TE_{10} H Field at the Port', 'FontSize', 20)
xlabel('x length, m', 'FontSize', 16)
ylabel('z length, m', 'FontSize', 16)

% TE 20

figure(7)
quiverC2Dfm(X,Y,real(H_x20),real(H_y20))
xlim([0, 0.075])
ylim([0, 0.0375])
title('Waveguide TE_{20} H Field at the Port', 'FontSize', 20)
xlabel('x length, m', 'FontSize', 16)
ylabel('z length, m', 'FontSize', 16)

% TE 01

figure(8)
quiverC2Dfm(X,Y,real(H_x01),real(H_y01))
xlim([0, 0.075])
ylim([0, 0.0375])
title('Waveguide TE_{01} H Field at the Port', 'FontSize', 20)
xlabel('x length, m', 'FontSize', 16)
ylabel('z length, m', 'FontSize', 16)

% TE 11

figure(9)
quiverC2Dfm(X,Y,real(H_x11),real(H_y11))
xlim([0, 0.075])
ylim([0, 0.0375])
title('Waveguide TE_{11} H Field at the Port', 'FontSize', 20)
xlabel('x length, m', 'FontSize', 16)
ylabel('z length, m', 'FontSize', 16)

% TE 21

figure(10)
quiverC2Dfm(X,Y,real(H_x21),real(H_y21))
xlim([0, 0.075])
ylim([0, 0.0375])
title('Waveguide TE_{21} H Field at the Port', 'FontSize', 20)
xlabel('x length, m', 'FontSize', 16)
ylabel('z length, m', 'FontSize', 16)



