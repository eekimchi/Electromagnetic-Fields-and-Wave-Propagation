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
z = linspace(0.113067/2, 0.113067/2, 30);

[X,Y] = meshgrid(x,y);

E_z = zeros(size(Y));

% TE101

E_x101 = cos(beta_x(2).*X).*sin(beta_y(1).*Y).*sin(beta_z(1).*z);

E_y101 = -sin(beta_x(2).*X).*cos(beta_y(1).*Y).*sin(beta_z(1).*z);

figure(1)
quiverC2D(X,Y,E_x101,E_y101)
xlim([0, 0.075])
ylim([0, 0.0375])
title('Resonant Cavity TE_{101} E Field at the XZ Plane', 'FontSize', 20)
xlabel('x length in m', 'FontSize', 16)
ylabel('z length in m', 'FontSize', 16)

% TE102

E_x102 = cos(beta_x(2).*X).*sin(beta_y(1).*Y).*sin(beta_z(2).*z);

E_y102 = -sin(beta_x(2).*X).*cos(beta_y(1).*Y).*sin(beta_z(2).*z);

figure(2)
quiverC2D(X,Y,E_x102,E_y102)
xlim([0, 0.075])
ylim([0, 0.0375])
title('Resonant Cavity TE_{102} E Field at the XZ Plane', 'FontSize', 20)
xlabel('x length in m', 'FontSize', 16)
ylabel('z length in m', 'FontSize', 16)

% TE201

E_x201 = cos(beta_x(3).*X).*sin(beta_y(1).*Y).*sin(beta_z(1).*z);

E_y201 = -sin(beta_x(3).*X).*cos(beta_y(1).*Y).*sin(beta_z(1).*z);

figure(3)
quiverC2D(X,Y,E_x201,E_y201)
xlim([0, 0.075])
ylim([0, 0.0375])
title('Resonant Cavity TE_{201} E Field at the XZ Plane', 'FontSize', 20)
xlabel('x length in m', 'FontSize', 16)
ylabel('z length in m', 'FontSize', 16)

% TE 011

E_x011 = cos(beta_x(1).*X).*sin(beta_y(2).*Y).*sin(beta_z(1).*z);

E_y011 = -sin(beta_x(1).*X).*cos(beta_y(2).*Y).*sin(beta_z(1).*z);

figure(4)
quiverC2D(X,Y,E_x011,E_y011)
xlim([0, 0.075])
ylim([0, 0.0375])
title('Resonant Cavity TE_{011} E Field at the XZ Plane', 'FontSize', 20)
xlabel('x length in m', 'FontSize', 16)
ylabel('z length in m', 'FontSize', 16)

% TE 103

E_x103 = cos(beta_x(2).*X).*sin(beta_y(1).*Y).*sin(beta_z(3).*z);

E_y103 = -sin(beta_x(2).*X).*cos(beta_y(1).*Y).*sin(beta_z(3).*z);

figure(5)
quiverC2D(X,Y,E_x103,E_y103)
xlim([0, 0.075])
ylim([0, 0.0375])
title('Resonant Cavity TE_{103} E Field at the XZ Plane', 'FontSize', 20)
xlabel('x length in m', 'FontSize', 16)
ylabel('z length in m', 'FontSize', 16)

% XZ Plane

x = linspace(0, 0.075, 30);
y = linspace(0.0375/2, 0.0375/2, 30);
z = linspace(0, 0.113067, 30);
z = z';

[X,Z] = meshgrid(x,z);

% TE 101

E_y101 = -(beta_x(2)./e0).*sin(beta_x(2).*X).*cos(beta_y(1).*y).*sin(beta_z(1).*Z);

figure(6)
contourf(X,Z,abs(E_y101))
colorbar 
colormap jet
xlim([0, 0.075])
ylim([0, 0.113067])
title('Resonant Cavity TE_{101} E_{Mag} Field at the XY Plane', 'FontSize', 20)
xlabel('x length in m', 'FontSize', 16)
ylabel('y length in m', 'FontSize', 16)

% TE 102

E_y102 = -(beta_x(2)./e0).*sin(beta_x(2).*X).*cos(beta_y(1).*y).*sin(beta_z(2).*Z);

figure(7)
contourf(X,Z,abs(E_y102))
colorbar 
colormap jet
xlim([0, 0.075])
ylim([0, 0.113067])
title('Resonant Cavity TE_{102} E_{Mag} Field at the XY Plane', 'FontSize', 20)
xlabel('x length in m', 'FontSize', 16)
ylabel('y length in m', 'FontSize', 16)

% TE 201

E_y201 = -(beta_x(3)./e0).*sin(beta_x(3).*X).*cos(beta_y(1).*y).*sin(beta_z(1).*Z);

figure(8)
contourf(X,Z,abs(E_y201))
colorbar 
colormap jet
xlim([0, 0.075])
ylim([0, 0.113067])
title('Resonant Cavity TE_{201} E_{Mag} Field at the XY Plane', 'FontSize', 20)
xlabel('x length in m', 'FontSize', 16)
ylabel('y length in m', 'FontSize', 16)

% TE 011

E_x011 = cos(beta_x(1).*X).*sin(beta_y(2).*y).*sin(beta_z(1).*Z);

figure(9)
quiverC2D(X,Z,E_x011,E_z)
xlim([0, 0.075])
ylim([0, 0.113067])
title('Resonant Cavity TE_{011} E Field at the XY Plane', 'FontSize', 20)
xlabel('x length in m', 'FontSize', 16)
ylabel('y length in m', 'FontSize', 16)

% TE 103

E_y103 = -(beta_x(2)./e0).*sin(beta_x(2).*X).*cos(beta_y(1).*y).*sin(beta_z(3).*Z);

figure(10)
contourf(X,Z,abs(E_y103))
colorbar 
colormap jet
xlim([0, 0.075])
ylim([0, 0.113067])
title('Resonant Cavity TE_{103} E_{Mag} Field at the XY Plane', 'FontSize', 20)
xlabel('x length in m', 'FontSize', 16)
ylabel('y length in m', 'FontSize', 16)


% ZY Plane

x = linspace(0.075/2, 0.075/2, 30);
y = linspace(0, 0.0375, 30);
y = y';
z = linspace(0, 0.113067, 30);

[Z,Y] = meshgrid(z,y);

% TE 101

E_y101 = -cos(beta_y(1).*Y).*sin(beta_z(1).*Z);

figure(11)
quiverC2D(Z,Y,E_z,E_y101)
xlim([0, 0.113067])
ylim([0, 0.0375])
title('Resonant Cavity TE_{101} E Field at the YZ Plane', 'FontSize', 20)
xlabel('y length in m', 'FontSize', 16)
ylabel('z length in m', 'FontSize', 16)


% TE 102

E_y102 = -sin(beta_x(2).*x).*cos(beta_y(1).*Y).*sin(beta_z(2).*Z);

figure(12)
quiverC2D(Z,Y,E_z,E_y102)
xlim([0, 0.113067])
ylim([0, 0.0375])
title('Resonant Cavity TE_{102} E Field at the YZ Plane', 'FontSize', 20)
xlabel('y length in m', 'FontSize', 16)
ylabel('z length in m', 'FontSize', 16)

% TE 201

E_y201 = -sin(beta_x(3).*x).*cos(beta_y(1).*Y).*sin(beta_z(1).*Z);

figure(13)
quiverC2D(Z,Y,E_z,E_y201)
xlim([0, 0.113067])
ylim([0, 0.0375])
title('Resonant Cavity TE_{201} E Field at the YZ Plane', 'FontSize', 20)
xlabel('y length in m', 'FontSize', 16)
ylabel('z length in m', 'FontSize', 16)

% TE 011

E_x011 = (beta_y(2)./e0).*cos(beta_x(1).*x).*sin(beta_y(2).*Y).*sin(beta_z(1).*Z);

figure(14)
contourf(Z,Y,abs(E_x011))
colorbar 
colormap jet
xlim([0, 0.113067])
ylim([0, 0.0375])
title('Resonant Cavity TE_{011} E_{Mag} Field at the YZ Plane', 'FontSize', 20)
xlabel('y length in m', 'FontSize', 16)
ylabel('z length in m', 'FontSize', 16)

% TE 103

E_y103 = -sin(beta_x(2).*x).*cos(beta_y(1).*Y).*sin(beta_z(3).*Z);

figure(15)
quiverC2D(Z,Y,E_z,E_y103)
xlim([0, 0.113067])
ylim([0, 0.0375])
title('Resonant Cavity TE_{103} E Field at the YZ Plane', 'FontSize', 20)
xlabel('y length in m', 'FontSize', 16)
ylabel('z length in m', 'FontSize', 16)




