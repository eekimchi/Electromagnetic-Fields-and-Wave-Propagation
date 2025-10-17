clc
clear
close all

% 2.4 GHz Resonator Cavity

% Width
a = 0.075;

% Height
b = 0.0375;

% Depth
c = 0.113067;

% permittivity
e0 = 8.854e-12;

% permeability
mu0 = (4*pi)*(10^(-7));

% frequency
f = 5.7e9;

% convert frequency to omega
omega_f = 2.*pi.*f;

% beta
beta = omega_f.*sqrt(mu0.*e0);

% Transverse Electric TE

for m = 1:1:5
for n = 1:1:5
for p = 1:1:5
    frTE(m,n,p) = ((3e8)/(2*pi)).*sqrt((((m-1).*pi)./a).^2 + (((n-1).*pi)./b).^2 + ((p.*pi)./c).^2);
end
end
end

% Transverse Magnetic TM

for m = 1:1:5
for n = 1:1:5
for p = 1:1:5
    frTM(m,n,p) = ((3e8)/(2*pi)).*sqrt((((m).*pi)./a).^2 + (((n).*pi)./b).^2 + (((p-1).*pi)./c).^2);
end
end
end

for m = 1:1:5
    beta_x(m) = ((m-1).*pi)./a;
end

for n = 1:1:5
    beta_y(n) = ((n-1).*pi)./b;
end

for p = 1:1:5
    beta_z(p) = (p.*pi)./c;
end

% XY Plane

x = linspace(0, 0.075, 30);
y = linspace(0, 0.0375, 30);
y = y';
z = 0.113067/2;

[X,Y] = meshgrid(x,y);

% Vectors
H_x103 = 1j.*sin(beta_x(2).*X).*cos(beta_y(1).*Y).*cos(beta_z(3).*z);
H_y103 = 1j.*cos(beta_x(2).*X).*sin(beta_y(1).*Y).*cos(beta_z(3).*z);

% Magnitude
H_x103m = 1j.*((beta_x(2).*beta_z(3))./(omega_f.*mu0.*e0)).*sin(beta_x(2).*X).*cos(beta_y(1).*Y).*cos(beta_z(3).*z);
H_z103 = -1j.*((-beta_z(3).^2 + beta.^2)./(omega_f.*mu0.*e0)).*cos(beta_x(2).*X).*cos(beta_y(1).*Y).*sin(beta_z(3).*z);

H_z103 = sqrt((H_x103m.^2)+(H_z103.^2));

figure(1)
contourf(X,Y,abs(H_z103),'edgecolor','none')
colorbar 
colormap jet
xlim([0, 0.075])
ylim([0, 0.0375])
title('Resonant Cavity TE_{103} H_{Mag} Field at the XZ Plane', 'FontSize', 20)
xlabel('x length, m', 'FontSize', 16)
ylabel('z length, m', 'FontSize', 16)
figure(2)
quiverC2Dm(X,Y,abs(H_x103),abs(H_y103))
xlim([0, 0.075])
ylim([0, 0.0375])
title('Resonant Cavity TE_{103} H_{Vector} Field at the XZ Plane', 'FontSize', 20)
xlabel('x length, m', 'FontSize', 16)
ylabel('z length, m', 'FontSize', 16)

% XZ Plane

x = linspace(0, 0.075, 30);
z = linspace(0, 0.113067, 30);
z = z';
y = 0.0375/2;

[X,Z] = meshgrid(x,z);

% Vectors
H_x103 = 1j.*sin(beta_x(2).*X).*cos(beta_y(1).*y).*cos(beta_z(3).*Z);
H_z103v = -1j.*cos(beta_x(2).*X).*cos(beta_y(1).*y).*sin(beta_z(3).*Z);

% Magnitude
H_x103m = 1j.*((beta_x(2).*beta_z(3))./(omega_f.*mu0.*e0)).*sin(beta_x(2).*X).*cos(beta_y(1).*y).*cos(beta_z(3).*Z);
H_z103 = -1j.*((-beta_z(3).^2 + beta.^2)./(omega_f.*mu0.*e0)).*cos(beta_x(2).*X).*cos(beta_y(1).*y).*sin(beta_z(3).*Z);

H_z103 = sqrt((H_x103m.^2)+(H_z103.^2));

figure(3)
contourf(X,Z,abs(H_z103),'edgecolor','none')
colorbar 
colormap jet
xlim([0, 0.075])
ylim([0, 0.113067])
title('Resonant Cavity TE_{103} H_{Mag} Field at the XY Plane', 'FontSize', 20)
xlabel('x length, m', 'FontSize', 16)
ylabel('y length, m', 'FontSize', 16)
figure(4)
quiverC2Dm(X,Z,imag(H_x103),imag(H_z103v))
xlim([0, 0.075])
ylim([0, 0.113067])
title('Resonant Cavity TE_{103} H_{Vector} Field at the XY Plane', 'FontSize', 20)
xlabel('x length, m', 'FontSize', 16)
ylabel('y length, m', 'FontSize', 16)

% YZ Plane

x = 0.075/2;
z = linspace(0, 0.113067, 30);
y = linspace(0, 0.0375, 30);
y = y';

[Z,Y] = meshgrid(z,y);

% Vectors
H_z103v = -1j.*cos(beta_x(2).*x).*cos(beta_y(1).*Y).*sin(beta_z(3).*Z);
H_y103 = 1j.*cos(beta_x(2).*x).*sin(beta_y(1).*Y).*cos(beta_z(3).*Z);

% Magnitude
H_x103m = 1j.*((beta_x(2).*beta_z(3))./(omega_f.*mu0.*e0)).*sin(beta_x(2).*x).*cos(beta_y(1).*Y).*cos(beta_z(3).*Z);
H_z103 = -1j.*((-beta_z(3).^2 + beta.^2)./(omega_f.*mu0.*e0)).*cos(beta_x(2).*x).*cos(beta_y(1).*Y).*sin(beta_z(3).*Z);

H_z103 = sqrt((H_x103m.^2)+(H_z103.^2));

figure(5)
contourf(Z,Y,abs(H_z103),'edgecolor','none')
colorbar 
colormap jet
xlim([0, 0.113067])
ylim([0, 0.0375])
title('Resonant Cavity TE_{103} H_{Mag} Field at the YZ Plane', 'FontSize', 20)
xlabel('y length, m', 'FontSize', 16)
ylabel('z length, m', 'FontSize', 16)
% figure(6)
% quiverC2Dm(Z,Y,imag(H_z103v),(H_y103))
% xlim([0, 0.113067])
% ylim([0, 0.0375])
% title('Resonant Cavity TE_{103} H_{Vector} Field at the YZ Plane', 'FontSize', 20)
% xlabel('x length, m', 'FontSize', 16)
% ylabel('z length, m', 'FontSize', 16)
