clc
clear
close all

% ========== CAVITY ==========

% Microwave dimensions mmmmmmmmmmmmmmmmm

c = 21*0.0254; % width
b = 9.5*0.0254; % height
a = 14*0.0254; % depth

% Transverse Electric TE

for m = 1:1:5
for n = 1:1:5
for p = 1:1:5
    microwavefrTE(m,n,p) = ((3e8)/(2*pi)).*sqrt((((m-1).*pi)./a).^2 + (((n-1).*pi)./b).^2 + ((p.*pi)./c).^2);
end
end
end

% Transverse Magnetic TM

for m = 1:1:5
for n = 1:1:5
for p = 1:1:5
    microwavefrTM(m,n,p) = ((3e8)/(2*pi)).*sqrt((((m).*pi)./a).^2 + (((n).*pi)./b).^2 + (((p-1).*pi)./c).^2);
end
end
end