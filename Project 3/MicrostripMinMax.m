clc
clear
close all


mu0 = 4*pi*10^(-7);
e0 = 8.85e-12;
er = 2.2;
epsilon = e0*er;

% 2.4 GHz

f = 2.4e9;

lambda = 1/(f*sqrt(mu0*epsilon));

open = 0.5 - 0.251;
opened = open*lambda

close = open + 0.25;
closed = close*lambda

% 4 GHz

f = 4e9;

lambda = 1/(f*sqrt(mu0*epsilon));

open = 0.5 - 0.251;
opened = open*lambda

close = open + 0.25;
closed = close*lambda


% 6 GHz

f = 6e9;

lambda = 1/(f*sqrt(mu0*epsilon));

open = 0.5 - 0.251;
opened = open*lambda

close = open + 0.25;
closed = close*lambda


