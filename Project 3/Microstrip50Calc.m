clc
clear
close all

Zo = 50;
er = 2.2;

% Pozar
% B = (377*pi)/(2*Zo*sqrt(er));
% 
% Wd = (2/pi)*(B - 1 - log(2*B - 1) + ((er - 1)/(2*er))*(log(B - 1) + 0.39 - ((0.61)/(er))));
% 
% W = Wd*1.575;
% 
% a = 1/(2*pi);
% 
% % this is to double check
% dW = 1/Wd;
% 
% ee = ((er + 1)/(2)) + ((er - 1)/(2))*((1)/(sqrt(1 + 12*dW)));
% 
% Zo = (120*pi)/(sqrt(ee)*(Wd + 1.393 + 0.667*log(Wd + 1.444)));
% 
% % A = ((Zo/60) * sqrt((er + 1)/ 2)) + ((er - 1)/(er + 1) * (0.23 + (0.11/er)));
% % 
% % Wd = (8*exp(A))/(exp(2*A) - 2);


% Balanis
w = 4.7e-3;
h = 1.575e-3;
t = 17.5e-6;

weffh = (w/h) + ((1.25/pi)*(t/h)*(1 + log((2*h)/t)));

hweff = 1/weffh;

ereff = ((er + 1)/2) + (((er-1)/2)*((1 + 12*weffh))^(-1/2));

Zc = ((120*pi)/sqrt(ereff))/(weffh + 1.393 + 0.667*log(weffh + 1.444));



